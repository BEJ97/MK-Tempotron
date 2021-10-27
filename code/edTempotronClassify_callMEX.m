function [weights, correctRate, CorrectFireRate, CorrectNonFireRate,epoch]= ...
    edTempotronClassify_callMEX(weights, IsTraining, SimWhichSet, datafolder, MAXEPO, lmd,targetRate,a,b)
% event-driven simulation of tempotron neurons,
% Zhao Bo, NTU, Singapore. 2013.
% Ref: Tempotron: a neuron that learns spike timing-based decisions
% Rober Gutig 2006 Nature Neuroscience

if ~exist('datafolder','var')
    datafolder = pwd;
end

if ~exist('MAXEPO','var'), MAXEPO = 1000; end
if (IsTraining), 
    maxEpoch = MAXEPO; 
    load ([datafolder,'/','PtnCellTrnTst_spk'], 'PtnCellTrn')
    PtnCell = PtnCellTrn;
    clear PtnCellTrn;
else
    maxEpoch = 1;   % if not training, only run 1 epoch
    if isequal(SimWhichSet,'training set')
        load ([datafolder,'/','PtnCellTrnTst_spk'], 'PtnCellTrn')
        PtnCell = PtnCellTrn;
        clear PtnCellTrn;
    elseif isequal(SimWhichSet,'testing set')
        load ([datafolder,'/','PtnCellTrnTst_spk'], 'PtnCellTst')
        PtnCell = PtnCellTst;
        clear PtnCellTst;
    else
        error('Error: SimWhichSet can only be ''training set'' or ''testing set''!');
    end
end

if ~exist('lmd','var')
    lmd = 2e-2;
end
if ~exist('targetRate','var')
    targetRate = 1;
end


use_single_exponential = 0 ;
show_figures = 0;
show_which_pattern = [1,1,1,1];  %[iImage, pp, neuron, indNeuronPerOutput];


V_thr = 1  ; V_rest = 0; %改
T = 1; %256;            % pattern duration
tau_m = 0.1; %20;         % tau_m
tau_s = tau_m/4;
dt = tau_m/1000;             % time tick for lookup table
mu = 0;  % momentum factor (used for accelerating learning); this factor is not used since mu=0

tau1=tau_m;  
tau2=tau_s;

% lookup table
t = 0:dt:T;
lut1 = exp(-t/tau1);
if ~use_single_exponential,
    lut2 = exp(-t/tau2);
end

% normalization coefficient.
% V0 = 1/max(exp(-(0:dt:5*tau1)/tau1)-exp(-(0:dt:5*tau1)/tau2)); 
if ~use_single_exponential,
    tmp = round(5*tau1/dt);
    V0 = 1/max( lut1(1:tmp) - lut2(1:tmp) ) +0;  %改1.8
else
    V0 = 1;
end



nOutputs = size(weights,2);
nNeuronPerOutput = size(weights,3);
nImages = length(PtnCell);
correctRate=zeros(1,maxEpoch);
dw_Past=zeros(size(weights));  % momentum for accelerating learning.


% -------- PtnCell to Arrays: (cell type is not supported in MATLAB coder)
[ SpkTimings,Addresses,TimeChnlLbl,NeuronsTgt, NeuronsOut, ...
    idx_convert, hasFieldTgt, hasFieldOut] = ...
    ConvertPtnCell2Arrays( PtnCell,nOutputs,nNeuronPerOutput );
% ----------------

mexBuiltOnTrnSet = 0;
mexBuiltOnTstSet = 0;

for epoch=1:maxEpoch
    numTotSlices = 0;
    numCorrectSlices = 0;
    numTotFireSlices = 0;
    numCorrectFireSlices = 0;
    numTotNonFireSlices = 0;
    numCorrectNonFireSlices = 0;
    
    
%         [NeuronsTgt,NeuronsOut,numCorrectSlices,numTotSlices,numCorrectFireSlices,...
%             numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,...
%             weights,dw_Past] = ...
%             edTempotronCore (nNeuronPerOutput,nOutputs,T,...
%             SpkTimings,Addresses,TimeChnlLbl,NeuronsTgt,NeuronsOut,...
%             use_single_exponential,tau1,V0,lut1,lut2,dt,V_thr,...
%             numCorrectSlices,numTotSlices,numCorrectFireSlices,...
%             numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,...
%             IsTraining,lmd,mu,weights,dw_Past);
    

    % ---------  speed up with mex --------------- 
    % We need two seperate mex files since train set and test set usually have differet sizes.
    % NOTE: IF YOU SOMEHOW FAIL TO BUILD MEX, COMMENT OUT THE TWO 'codegen'
    % LINES BELOW AND USE THE ORIGINAL MATLAB FUNCTION 'edTempotronCore' TO
    % REPLACE 'edTempotronCore_onTrnSet_mex' AND 'edTempotronCore_onTstSet_mex'.
    % IT RUNS MUCH SLOWER WITHOUT MEX!
    if isequal(SimWhichSet,'training set')
        % build mex:
        if mexBuiltOnTrnSet ==0
            fprintf('Building mex function: edTempotronCore_onTrnSet_mex ...\n');
            codegen edTempotronCore -args {nNeuronPerOutput,nOutputs,T,SpkTimings,Addresses,TimeChnlLbl,NeuronsTgt,NeuronsOut,use_single_exponential,tau1,V0,lut1,lut2,dt,V_thr,numCorrectSlices,numTotSlices,numCorrectFireSlices,numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,IsTraining,lmd,mu,weights,dw_Past,a,b} -o edTempotronCore_onTrnSet_mex
            mexBuiltOnTrnSet = 1;
        end
        % call mex function:
        [NeuronsTgt,NeuronsOut,numCorrectSlices,numTotSlices,numCorrectFireSlices,...
            numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,...
            weights,dw_Past] = ...
            edTempotronCore_onTrnSet_mex (nNeuronPerOutput,nOutputs,T,...
            SpkTimings,Addresses,TimeChnlLbl,NeuronsTgt,NeuronsOut,...
            use_single_exponential,tau1,V0,lut1,lut2,dt,V_thr,...
            numCorrectSlices,numTotSlices,numCorrectFireSlices,...
            numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,...
            IsTraining,lmd,mu,weights,dw_Past,a,b);
    elseif isequal(SimWhichSet,'testing set')
        % build mex:
        if mexBuiltOnTstSet==0
            fprintf('Building mex function: edTempotronCore_onTstSet_mex ...\n');
            codegen edTempotronCore -args {nNeuronPerOutput,nOutputs,T,SpkTimings,Addresses,TimeChnlLbl,NeuronsTgt,NeuronsOut,use_single_exponential,tau1,V0,lut1,lut2,dt,V_thr,numCorrectSlices,numTotSlices,numCorrectFireSlices,numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,IsTraining,lmd,mu,weights,dw_Past,a,b} -o edTempotronCore_onTstSet_mex
            mexBuiltOnTstSet = 1;
        end
        % call mex function:
        [NeuronsTgt,NeuronsOut,numCorrectSlices,numTotSlices,numCorrectFireSlices,...
            numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,...
            weights,dw_Past] = ...
            edTempotronCore_onTstSet_mex (nNeuronPerOutput,nOutputs,T,...
            SpkTimings,Addresses,TimeChnlLbl,NeuronsTgt,NeuronsOut,...
            use_single_exponential,tau1,V0,lut1,lut2,dt,V_thr,...
            numCorrectSlices,numTotSlices,numCorrectFireSlices,...
            numTotFireSlices,numTotNonFireSlices,numCorrectNonFireSlices,...
            IsTraining,lmd,mu,weights,dw_Past,a,b);
    end

    

    % weights after each epoch's modification
    TrainedWt=weights;
    correctRate(epoch)= numCorrectSlices / (numTotSlices);
    
    CorrectFireRate = numCorrectFireSlices/numTotFireSlices;
    CorrectNonFireRate = numCorrectNonFireSlices/numTotNonFireSlices;

    % Actually, correctRate(epoch) is not the exact rate (except when it reaches 1), 
    % it's just a rough estimation since the weight is modified for each pattern during this epoch.
    
    if(IsTraining), save([datafolder,'/','TrainedWt'],'TrainedWt','correctRate','epoch'); end
    if(IsTraining), timedLog( ['epoch: ', num2str(epoch), ', correct Fire rate: ', num2str(CorrectFireRate),...
            ', correct NonFire rate: ', num2str(CorrectNonFireRate), ', total correct rate: ', num2str(correctRate(epoch))] ); end
    
    if(IsTraining),
        Rates = correctRate(max(1,epoch-9):epoch);
        avgRate(epoch) = mean(Rates);
        [tmp idx1] = sort( -1* avgRate(max(1,epoch-9):epoch) );
        condition1 = correctRate(epoch)==1;         % all correct
        condition2 = sum(Rates==Rates(1)) ==10;     % no change of correctRate for 10 consecutive epochs
        condition3 = isequal(idx1,1:10);            % avgRate decreases for 10 consecutive epochs
        condition4 = (epoch>10) & ( sum(abs(Rates-mean(Rates))) < 1e-3 ); % rate almost saturates
        condition5 = correctRate(epoch)>= targetRate;
        if condition1 || condition5 
            if condition1, timedLog('Training ends: 100%% rate achieved'); end
            if condition5, timedLog(sprintf('target rate %.6f achieved',targetRate)); end
            correctRate = correctRate(1:epoch);
            break;     % break, no need to run more epochs.
        end
		if condition2 || condition3 || condition4
			% reduce the learning rate.
            if lmd> 1e-6
                lmd = lmd/2;
                timedLog( sprintf('lmd changed to %.6f',lmd) );
            else
                timedLog( 'corr rate almost saturates. stop training');
                correctRate = correctRate(1:epoch);
                break
            end
		end
        if epoch==maxEpoch, timedLog(['Training ends: maxEpoch (' num2str(maxEpoch) ') achieved']); end
    end
   
    
end	% end of one epoch




% ---------- Arrays to PtnCell :
hasFieldTgt=1;
hasFieldOut=1;
[ PtnCell ] = ConvertArrays2PtnCell(SpkTimings,Addresses,...
    TimeChnlLbl,NeuronsTgt,NeuronsOut,idx_convert,hasFieldTgt,hasFieldOut);
% ----------------




if ~IsTraining, 
    if isequal(SimWhichSet,'training set')
        PtnCellTrn = PtnCell;
        clear PtnCell;
        save ([datafolder,'/','PtnCellTrn_out'], 'PtnCellTrn')
    elseif isequal(SimWhichSet,'testing set')
        PtnCellTst = PtnCell;
        clear PtnCell;
        save ([datafolder,'/','PtnCellTst_out'], 'PtnCellTst')
    end
end








