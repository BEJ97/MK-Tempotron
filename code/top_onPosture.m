
% clear;clc;
close all
format compact
zonghe = zeros(10,1);
for i = 1:1
reRunAll = 0;   % if 1, rerun all steps; if 0, some steps may be skipped if corresponding mat files exist.
reGenCINs     = 0;
reGetTHs      = 0;
reGenPtnCell  = 0;
reSplitDataset= 0;
reSpkConvert  = 0;
reInitWeights = 0;
reTrnWeights  = 0;
continueTrnWts= 0;
reSimulation  = 1;
reAnalyzeResults = 1 ;

a = 0.5;
b = 0.5;

maxEpoch = 15; 

lmd = 1e-1;
% lmd = 2e-2;
% lmd = 1e-2;
% lmd = 1e-3;
% lmd = 1e-4;
% lmd = 1e-5;

nNeuronPerOutput = 10;      % number of neurons at each output channel. (redundancy, majority voting)

nGroup = 3;
size1 = 32;
size2 = 32;

% integ switch
tau_m = 20e6;
timeslice = 1.5*tau_m;    	% for reconstructing img from event slice
% conv_disp_interval = 10*timeslice; 
conv_disp_interval = timeslice; 

% filters:
minFS = 3; maxFS = 9; sSFS  = 2;
rot = [0 45 90 135];
circularRField   = 1;

useSmooth = 1;
if useSmooth==1, str1='SM'; else str1=''; end


CONV_refresh_window_only = 1;  % if 1, only refresh a window. otherwise, refresh whole S1 map.

% Leakage rate (on the fly conv)
cteloss = 5e-8;     % 1/(20e3),  leakage per nanosecond.   

% MAX
doMAXforEachEvent = 0;      
% if 1, MAX operation is performed for each event. (only for the pixels
% within the region affected by convolution.); if 0, MAX is performed only 
% after peak detected. (for whole map).
if doMAXforEachEvent==1, maxStr='_myWinMAX'; else maxStr='_myMAX'; end

post_proc_th     = 0.2;%0.8; %0.5; 

% whether or not to use centroidShift on C1 responses.
useCentroidShift = 1;

showOnTheFlyConv = 0;
showLayerFigures = 0;

% ratio of training samples to total
ratio_Trn = 0.8;

% data folder:
folder0 = ['dnew_posture_sz',num2str(size1), '_4x4filt',str1,maxStr,num2str(post_proc_th),...
    '_wRedun',num2str(nNeuronPerOutput) ,'_tau',mynum2str(tau_m,'%0.0e'),'_cte',mynum2str(cteloss)];
if ~exist(folder0,'dir'), mkdir (folder0); end



%% prepare lut for exp(-t/tau), used for leaky integration cell.

use_single_exponential = 0;   %改

dt = tau_m/1e3; % dt for lut
tau_s = tau_m/4;
tau1 = tau_m;
tau2 = tau_s;

% lookup table
t = 0:dt:10*tau1;
lut1 = exp(-t/tau1);
lut2 = exp(-t/tau2);

% normalization coefficient.
% V0 = 1/max(exp(-(0:dt:5*tau1)/tau1)-exp(-(0:dt:5*tau1)/tau2)); 
if ~use_single_exponential,
    tmp = round(5*tau1/dt);
    V0 = 1/max( lut1(1:tmp) - lut2(1:tmp) );
else
    V0 = 1;
end




%%

CINsfile = [folder0, '/', 'posture.mat'];
if ( (reRunAll) || (reGenCINs) || ~exist(CINsfile,'file') )
    [NumCINs]=prepareCINs(folder0);
end

load (CINsfile,'NumCINs')   % note: CIN in ns





%%

file = [folder0, '/', 'THs.mat'];
if ( (reRunAll) || (reGetTHs) || ~exist(file,'file') )
    [ THs ] = getThStats(CINsfile, V0, dt,use_single_exponential, lut1,lut2, nGroup, folder0);
end
load (file)



file = [folder0, '/', 'PtnCell500.mat'];
if ( (reRunAll) || (reGenPtnCell) || ~exist(file,'file') )
    timedLog('generating PtnCell...');
    
    PtnCell = cell(1, NumCINs);        % data to tempotron. one cell for each CIN.
    
    load (CINsfile,'CINs','Labels');
    
    parfor ind = 1:NumCINs

        disp(ind);

        CIN = CINs{ind};
        label = double(Labels(ind));

        % ------------------------------------------------------------------------
        % motion symbol detector (one leaky integration neuron and peak detection)
        [ Timings, V, t ] = MotionSymbolDetector_withCentroid( CIN, ...
            V0,dt,use_single_exponential,lut1,lut2, THs, label, nGroup, size1, 0) ;
        
        % -----------------------------------------------------------------------

        % params initialization
        [ paramCell,stateCell,sqfilters,filtersize,numScales,numRot,numFilters] = ...
            init_paramCell_stateCell( minFS,maxFS,sSFS,rot,circularRField,useSmooth, cteloss,size1,size2 );

        % on the fly convolution and MAX operation
        [ AllVec,Time_Chnl_Lbl ] = conv_and_max( CIN,numScales,numRot,numFilters,size1,size2,...
            Timings,paramCell,stateCell, timeslice,...
            sqfilters,filtersize, doMAXforEachEvent, post_proc_th,...
            showOnTheFlyConv,showLayerFigures,rot,circularRField,conv_disp_interval,...
            CONV_refresh_window_only, useCentroidShift) ;

	
		[AllVec,AllAddr] = sort(AllVec,1,'descend');      
        
        % -------------
        nkeep = 100;
        
        PtnCell{ind}.AllVec  = AllVec(1:nkeep,:);
        PtnCell{ind}.AllAddr = AllAddr(1:nkeep,:);
        
        PtnCell{ind}.Time_Chnl_Lbl = Time_Chnl_Lbl;

        %toc
    end

    save (file, 'PtnCell')
end

clear CINs 
clear Images Labels
clear PtnCell

load (file)

% maxResp = findMAXofAllPtn(PtnCell);



%%

% --- random split
file = [folder0, '/', 'PtnCellTrnTst_raw.mat'];
if ( (reRunAll) || (reSplitDataset) || ~exist(file,'file') )
    [PtnCellTrn, PtnCellTst, indTrn, indTst, indNaN] = RandSplit(PtnCell,ratio_Trn,nGroup);
    save (file, 'PtnCellTrn', 'PtnCellTst', 'indTrn', 'indTst', 'indNaN')
end
load (file)
maxResp = findMAXofAllPtn(PtnCellTrn);

clear PtnCell
clear PtnCellTrn PtnCellTst indTrn indTst indNaN



% --- C1 resp to spike conversion
file = [folder0, '/', 'PtnCellTrnTst_spk.mat'];
if ( (reRunAll) || (reSpkConvert) || ~exist(file,'file') )
    myC1respToSpk(maxResp,folder0);
end





%% --- Tempotron
% initialize weights
nAfferents = size1*size2*(length(minFS:sSFS:maxFS)*length(rot));  
nOutputs = nGroup;
mu_init_wt = 0;
sigma_init_wt = 0.1;   % empirical 4/nAfferents=0.005
file = [folder0,'/','weights0.mat'];
if ( (reRunAll) || (reInitWeights) || ~exist(file,'file') )
    timedLog('initiating weights...');
    weights = mu_init_wt + sigma_init_wt * randn(nAfferents, nOutputs, nNeuronPerOutput);       % with redundancy
    save (file, 'weights')
end
load (file)

% TRAINING
tic
IsTraining = 1;
% maxEpoch = 10; 
SimWhichSet = 'training set';
file = [folder0,'/','TrainedWt.mat'];
TrnWtsFromInit = ( (reRunAll) || (reTrnWeights) || ~exist(file,'file') );
if ( (TrnWtsFromInit) || (continueTrnWts) )
    if ( TrnWtsFromInit )
        timedLog('start training ..');
 %        [TrainedWt,correctRate, tmp2, tmp3,epoch] = edTempotronClassify(weights, ...
 %            IsTraining, SimWhichSet, folder0, maxEpoch, lmd, targetRate,a,b);
        [TrainedWt,correctRate, tmp2, tmp3,epoch] = edTempotronClassify_callMEX(weights, ...
            IsTraining, SimWhichSet, folder0, maxEpoch, lmd,1,a,b);
    elseif (continueTrnWts)
        timedLog('continue training ..');
        load (file)
%         [TrainedWt,correctRate, tmp2, tmp3,epoch] = edTempotronClassify(TrainedWt, ...
%             IsTraining, SimWhichSet, folder0,maxEpoch, lmd,targetRate,a,b);
        [TrainedWt,correctRate, tmp2, tmp3,epoch] = edTempotronClassify_callMEX(TrainedWt, ...
            IsTraining, SimWhichSet, folder0,maxEpoch, lmd,1,a,b);
    end
    timeTrain = toc/60; % min
    if timeTrain<60
        timedLog(['Training finished, time taken: ',num2str(timeTrain),' min'])
    else
        timedLog(['Training finished, time taken: ',num2str(timeTrain/60), ' hrs'])
    end
    save (file, 'TrainedWt','correctRate','epoch')
end
load (file)



% SIMULATION
tic
IsTraining = 0;
file = [folder0,'/','RawResults.mat'];
if ( (reRunAll) || (reSimulation) || ~exist(file,'file') || ...
        ~exist([folder0,'/','PtnCellTrn_out.mat'],'file') || ~exist([folder0,'/','PtnCellTst_out.mat'],'file') )
    timedLog('start simulation ..');
%     [TrainedWt, RR_trn, RR_trn_TP, RR_trn_TN] = edTempotronClassify(TrainedWt, ...
%         IsTraining, 'training set',folder0, 1, lmd);
[TrainedWt, RR_trn, RR_trn_TP, RR_trn_TN] = edTempotronClassify_callMEX(TrainedWt, ...
        IsTraining, 'training set',folder0, 1, lmd,1,a,b);
%     [TrainedWt, RR_tst, RR_tst_TP, RR_tst_TN] = edTempotronClassify(TrainedWt, ...
%         IsTraining, 'testing set',folder0, 1, lmd,1,a,b);
    
    [TrainedWt, RR_tst, RR_tst_TP, RR_tst_TN] = edTempotronClassify_callMEX(TrainedWt, ...
        IsTraining, 'testing set',folder0, 1, lmd,1,a,b);
    timeSim = toc/60; % min
    if timeSim<60
        timedLog(['Sim finished, time taken: ',num2str(timeSim),' min'])
    else
        timedLog(['Sim finished, time taken: ',num2str(timeSim/60), ' hrs'])
    end

    save (file, 'RR_trn', 'RR_trn_TP', 'RR_trn_TN', 'RR_tst', 'RR_tst_TP', 'RR_tst_TN');
end
load (file)
fprintf('======================================================\n');
fprintf('Raw correct Rate (i.e. rate of correctly responded neurons for all slices):\n');
fprintf('\t training set: total succ rate %.2f %%, \t true positive %.2f %%, \t true negative %.2f %% \n',  ...
    RR_trn*100, RR_trn_TP*100, RR_trn_TN*100);
fprintf('\t testing set : total succ rate %.2f %%, \t true positive %.2f %%, \t true negative %.2f %% \n',  ...
    RR_tst*100, RR_tst_TP*100, RR_tst_TN*100);


% analyze raw results, use majority voting to make final decision.
file = [folder0,'/','FinalResults.mat'];
if ( (reRunAll) || (reAnalyzeResults) || ~exist(file,'file') )
    [ FinalCorrRate_Trn, FinalCorrRate_Tst, decision_Trn, decision_Tst, freqsCell_Trn, freqsCell_Tst ] = ...
        AnalyzeResults_MajorityVoting(nOutputs, nNeuronPerOutput, folder0);
    save (file, 'FinalCorrRate_Trn', 'FinalCorrRate_Tst', 'decision_Trn', ...
        'decision_Tst', 'freqsCell_Trn', 'freqsCell_Tst')
end
load (file)
fprintf('======================================================\n');
fprintf('Final correct Rate (i.e. succ rate of categorization, by majority voting):\n');
fprintf('\t training set:  %.2f %% \n',  FinalCorrRate_Trn*100);
fprintf('\t testing set :  %.2f %% \n',  FinalCorrRate_Tst*100);
fprintf('\n\n');
zonghe(i,1) = FinalCorrRate_Tst*100;
end




