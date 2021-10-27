function [weights, correctRate, CorrectFireRate, CorrectNonFireRate,epoch]= ...
    edTempotronClassify(weights, IsTraining, SimWhichSet, datafolder, MAXEPO, lmd,targetRate) 
%IsTraining=1 ，SimWhichSet = 'training set'， weights=3D矩阵 发放数*4(输出类别数)*10
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
        PtnCell = PtnCellTrn;   %训练集cell
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
    lmd = 2e-2;  %学习率
end
if ~exist('targetRate','var')
    targetRate = 1;
end


use_single_exponential = 0 ;
show_figures = 0;
show_which_pattern = [1,1,1,1];  %[iImage, pp, neuron, indNeuronPerOutput];


V_thr = 1; V_rest = 0;
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
    V0 = 1/max( lut1(1:tmp) - lut2(1:tmp) );
else
    V0 = 1;
end



nOutputs = size(weights,2);   % =4输出类别数
nNeuronPerOutput = size(weights,3);  % =10
nImages = length(PtnCell);  %样本数量
correctRate=zeros(1,maxEpoch);
dw_Past=zeros(size(weights));  % momentum for accelerating learning.


for epoch=1:maxEpoch
    numTotSlices = 0;
    numCorrectSlices = 0;
    numTotFireSlices = 0;
    numCorrectFireSlices = 0;
    numTotNonFireSlices = 0;
    numCorrectNonFireSlices = 0;
    
    %for iImage = 1:nImages
    for iImage = randperm(nImages)
    
        nSlices = size(PtnCell{iImage}.AllVec,2);      %AllVec列数即为切片数
        PtnCell{iImage}.Tgt = false(nNeuronPerOutput*nOutputs,nSlices);
        PtnCell{iImage}.Out = false(nNeuronPerOutput*nOutputs,nSlices);
        
        %for pp=1:nSlices
        for pp = randperm(nSlices)
            
            addr= PtnCell{iImage}.AllAddr(:,pp);			% spike address
            
            ptn = PtnCell{iImage}.AllVec(:,pp);				% spike timings
            
            % spike timing = T means original C1 =0. remove them if any
            ptn = ptn(ptn<T); 
            addr = addr(ptn<T);
            
            lbl = PtnCell{iImage}.Time_Chnl_Lbl(3,pp);
            tgt = false(nNeuronPerOutput*nOutputs,1);    % target
            tgt(lbl*nNeuronPerOutput+1: (lbl+1)*nNeuronPerOutput) = true;
            PtnCell{iImage}.Tgt(:,pp) = tgt;
			
			
			nAfferents = length(ptn);	% number of effective afferent. (e.g. 100)
			
            P = [ptn, addr, (1:nAfferents)'];

            if ~use_single_exponential,
                % peak is delayed by 0.462* tau_m. (FOR tau2=tau1/4.)
                % V0*(exp(-0.462*tau1/tau1)-exp(-0.462*tau1/tau2))
                peak_delay = 0.462*tau1;
                Pd = [min( P(:,1)+peak_delay, T), -1*ones(size(P(:,2))), -1*ones(size(P(:,2)))]; % peak timings.
                P = [P;Pd];
            end

            [~, idx_tmp] = sort(P(:,1),'ascend');
            P = P(idx_tmp,:,:);

            P = [P; T,-1,-1];    % the last dummy event does not require checking of firing.
			
			numEvt = size(P,1);
			
            
            for neuron = 1:nOutputs

                for indNeuronPerOutput = 1:nNeuronPerOutput
                    out = false;        % output

					% for monitoring the maximal value
                    Vmax= -inf; 
                    tmax= -inf;
                    t_fire = -inf;
					fired=false;   
                    if isequal(show_which_pattern, [iImage, pp, neuron, indNeuronPerOutput]) && (show_figures)
						Vm1=zeros(1,numEvt);  			% trace a pattern
                    end
					t_last = -1;
                    t_latestRealEvt = -inf;
					cnt_evts_of_same_timestamp = 1;
                    Vm = 0;
                    Vm_K1 =0;
                    Vm_K2 = 0;
					
					wts_eff =weights(addr,neuron,indNeuronPerOutput);	% effective weights.
					
					for i = 1:numEvt    
						
						t = P(i,1);
						addr_i = P(i,2);
						c = P(i,3);
						
						delta_t = t-t_last;

						% Different input spikes may have the same time of occurrence. Only
						% the last spike at that time triggers the possible output firing
						% and learning. As in the code, we check the firing and do the
						% learning for the last Vm when current delta_t is not 0. So there
						% is a delay of one event (to be exact, a delay of 1 to N events, N 
						% is the number of events that have the same timestamp at last tick).
                        % Doing this is only to make the event-driven output to be as
                        % similar as possible to the time-driven results.

						condition_lastVm_checkup = ((delta_t>0)&&(i>1)) || (i==numEvt);
                        if ~condition_lastVm_checkup,
							if i~=1
								cnt_evts_of_same_timestamp = cnt_evts_of_same_timestamp +1;
							end
                        else

							if Vm>Vmax
								Vmax=Vm; 
								tmax=t_last;
							end

							if Vm>=V_thr && fired==false % fire
								fired=true;
								t_fire=t_last;
								out =true;
								if use_single_exponential
									Vm = 1.2*Vm; % to make the output spike more noticeable
								end
							end

							if isequal(show_which_pattern, [iImage, pp, neuron, indNeuronPerOutput]) && (show_figures)
								for j= 1:cnt_evts_of_same_timestamp
									Vm1(i-j)=Vm;        % record Vm for last time tick (several events happened at that tick).
								end
							end

							cnt_evts_of_same_timestamp = 1;     % important to reset this counter
                        end
                        
                        refractory = fired ;%&& (t-t_fire <= t_ref);
                        if refractory,
                            break;
                        else

                            lut_addr = round(delta_t/dt)+1;
                            if lut_addr<=length(lut1)
                                Sc1 = lut1(lut_addr);
                            else
                                Sc1 = 0;
                            end
                            Vm_K1 = Sc1*Vm_K1;
                            if (c ~= -1)    % only for the original input spikes, not the dummy ones.
                                Vm_K1 = Vm_K1 + V0*weights(addr_i,neuron,indNeuronPerOutput);
                            end
                            if ~use_single_exponential,
                                if lut_addr<=length(lut2)
                                    Sc2 = lut2(lut_addr);
                                else
                                    Sc2 = 0;
                                end
                                Vm_K2 = Sc2*Vm_K2;
                                if (c ~= -1)
                                    Vm_K2 = Vm_K2 + V0*weights(addr_i,neuron,indNeuronPerOutput);
                                end
                                Vm = Vm_K1-Vm_K2;
                            else
                                Vm = Vm_K1; 
                            end

                            if c~=-1  % real event. not dummy ones.
                                t_latestRealEvt = t;
                            end

                            t_last = t;      
                        end
                        
					end % end of events

                    if Vmax<=0
                        tmax= t_latestRealEvt;
                    end
					
                    if isequal(show_which_pattern, [iImage, pp, neuron, indNeuronPerOutput]) && (show_figures)
						x = P(:,1)';
						figure(1); plot(x,Vm1,'-b',x, V_thr*ones(size(x)),'--r'); xlim([0,T]);
						title( strcat('stream#: ',num2str(iImage), '; slice#: ',num2str(pp), '; population#: ',num2str(neuron),'; neuron#: ',num2str(indNeuronPerOutput)) ); drawnow;
                    end
					
            
                    indTotNeuronOutput = ((neuron-1)*nNeuronPerOutput+indNeuronPerOutput);

                    PtnCell{iImage}.Out(indTotNeuronOutput, pp) = out;
                    numCorrectSlices = numCorrectSlices + double((tgt(indTotNeuronOutput)==out)) ;
                    numTotSlices = numTotSlices +1;

                    numCorrectFireSlices = numCorrectFireSlices +  double((tgt(indTotNeuronOutput)==1)&(tgt(indTotNeuronOutput)==out));
                    numTotFireSlices = numTotFireSlices + double((tgt(indTotNeuronOutput)==1));
                    numTotNonFireSlices = numTotNonFireSlices + double((tgt(indTotNeuronOutput)==0));
                    numCorrectNonFireSlices = numCorrectNonFireSlices + double((tgt(indTotNeuronOutput)==0)&(tgt(indTotNeuronOutput)==out));

					

                    % --------- TRAINING start ---------------
					if (IsTraining)  
						if out ~= tgt(indTotNeuronOutput)   % error (update weights) 
							
                            % calculate K_tmax, PSP kernel values of all channels at tmax.
                            index = ( P(:,1)<=tmax ) & ( P(:,2)~=-1 );
                            P1 = P(index,:);
                            % may have repeated addresses. i.e. multispike per afferent.

                            K_tmax = zeros(nAfferents,1);
                            K_tmax1 = zeros(nAfferents,1);
                            K_tmax2 = zeros(nAfferents,1);
                            
                            ts = P1(:,1);
                            addrs = P1(:,2);
                            cs = P1(:,3);
                            delta_ts = tmax - ts;
                            lut_addrs = round(delta_ts /dt)+1;
                            Sc1s = zeros(size(lut_addrs));
                            Sc1s(lut_addrs<=length(lut1)) = lut1(lut_addrs);
                            K_tmax1(cs) = V0*Sc1s;
                            if ~use_single_exponential,
                                Sc2s = zeros(size(lut_addrs));
                                Sc2s(lut_addrs<=length(lut2)) = lut2(lut_addrs);
                                K_tmax2(cs) = V0*Sc2s;
                                K_tmax = K_tmax1-K_tmax2;
                            else
                                K_tmax = K_tmax1; 
                            end
   
                            
							if fired==false    % LTP
								Dw=lmd*K_tmax;
							else           % LTD
								Dw=-1*lmd*K_tmax;
							end
							A1= wts_eff; 	%weights(addr,neuron,indNeuronPerOutput);
							dwPst=dw_Past(addr,neuron,indNeuronPerOutput);
							A1 = A1 + Dw + mu*dwPst.*(Dw~=0);
							weights(addr,neuron,indNeuronPerOutput)=A1; 	% ****
							dwPst(Dw~=0)=Dw(Dw~=0);
							dw_Past(addr,neuron,indNeuronPerOutput) = dwPst;
						end    
					end
                    % -------- TRAINING end ------------------



                end % end of one NeuronPerOutput
            
            end % end of one neuron (population)

        end % end of one slice
        
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








