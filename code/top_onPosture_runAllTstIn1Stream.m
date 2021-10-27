% run the algo on one piece of events. AllTstIn1Stream

top_onPosture



%% AllTstIn1CIN
folder1 = [folder0,'/','AllTstIn1'];
if ~exist(folder1, 'dir'), mkdir(folder1); end


file = [folder1,'/','CIN_AllTstIn1.mat'];

if ( (reRunAll) || (reSplitDataset) || ~exist(file,'file') )
    timedLog('combine all test CIN into 1 ...');

    load ([folder0, '/', 'posture.mat'], 'CINs','Labels');
    load ([folder0, '/', 'PtnCellTrnTst_raw.mat'], 'indTst')
    
    CINs = CINs(indTst);
    Labels = Labels(indTst);
    
    CIN_AllTstIn1 = [];
    tLabel_AllTstIn1 = [];
    % tstart0, label0
    % tend0,   label0
    % tstart1, label1
    % tend1,   label1 ...
    
    gapBTWN2 = 2.5*tau_m;       
    gapEmpty = 2.5*tau_m/10;    
  
    nRound = 2; % two rounds means 1 2 3, 1 2 3.
    
    inds = cell(nGroup,nRound);
    for i = 1:nGroup
        e = find(Labels == i-1);
        N = length(e);
        nEachRound = round(N/nRound);
        for iRound = 1:nRound
            if iRound~=nRound
                inds{i,iRound} = e( (iRound-1)*nEachRound+1 : iRound*nEachRound );
            else
                inds{i,iRound} = e( (iRound-1)*nEachRound+1 : N );
            end
        end
    end
    
    
    for iRound = 1:nRound
        for i = 1:nGroup
            %e = find(Labels == i-1);
            %N = length(e);
            %for j = 1:N

            e = inds{i,iRound} ;
            N = length(e);
            for j = 1:N

                CIN = CINs{e(j)};
                label = double( Labels(e(j)) );
                CIN(:,1) = CIN(:,1)-CIN(1,1);

                if isempty(CIN_AllTstIn1)
                    CIN_AllTstIn1 = [CIN_AllTstIn1; CIN];
                else
                    CIN(:,1) = CIN(:,1)+gapBTWN2+CIN_AllTstIn1(end,1);
                    numnoise = 20;
                    t2 = CIN(1,1);
                    t1 = CIN_AllTstIn1(end,1)+gapEmpty;
                    ts_noise = round( ( rand(numnoise,1)*(t2-t1)+t1 ) /100)*100;
                    ts_noise = sort(ts_noise);
                    x_noise = randi(size1,numnoise,1);
                    y_noise = randi(size2,numnoise,1);
                    CIN_noise = [ts_noise, zeros(numnoise,1), -1*ones(numnoise,1), x_noise,y_noise,ones(numnoise,1)];
                    CIN_AllTstIn1 = [CIN_AllTstIn1; CIN_noise; CIN];
                end

                if j==1
                    tmp1 = [CIN(1,1), label];
                    tLabel_AllTstIn1 = [tLabel_AllTstIn1; tmp1];
                elseif j==N
                    tmp1 = [CIN(end,1), label];
                    tLabel_AllTstIn1 = [tLabel_AllTstIn1; tmp1];
                end

            end
        end
    end
    
    clear CINs Labels
    
    save (file, 'CIN_AllTstIn1','tLabel_AllTstIn1')

end

load (file)

CIN = CIN_AllTstIn1;
tlabel = tLabel_AllTstIn1;

clear CIN_AllTstIn1 tLabel_AllTstIn1

% scatter3(CIN(:,1),CIN(:,4),CIN(:,5),'.'); 
% % ylim([1,size1]); zlim([1,size2]);
% xlabel('time');ylabel('x addr'), zlabel('y addr');






%%

file = [folder1,'/','Timings.mat'];

if ( (reRunAll) || (reSplitDataset) || ~exist(file,'file') )
    
    timedLog('integ AllTstIn1 ...');
    % ------------------------------------------------------------------------
    % motion symbol detector (one leaky integration neuron and peak detection)
    [ Timings, V, t ] = MotionSymbolDetector_withCentroid_AllTstIn1( CIN, ...
        V0,dt,use_single_exponential,lut1,lut2, THs, tlabel, nGroup, size1, 0) ;
    save (file, 'Timings')
end
load (file)


file = [folder1,'/','PtnCellTrnTst_raw.mat'];

if ( (reRunAll) || (reSplitDataset) || ~exist(file,'file') )
    timedLog('conv and max ...');
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

    numSurvived = sum(AllVec>0);

    ind = 1;
    nkeep = 100;
    PtnCell{ind}.AllVec  = AllVec(1:nkeep,:);
    PtnCell{ind}.AllAddr = AllAddr(1:nkeep,:);

    PtnCell{ind}.Time_Chnl_Lbl = Time_Chnl_Lbl;

    PtnCellTst = PtnCell;

    save (file, 'PtnCellTst')
    % clear PtnCell
    clear AllVec AllAddr Time_Chnl_Lbl
    clear stateCell sqfilters paramCell THs

end
load (file)



%%
% --- C1 resp to spike conversion
file = [folder1, '/', 'PtnCellTrnTst_spk.mat'];
if ( (reRunAll) || (reSplitDataset) || ~exist(file,'file') )
    timedLog('resp to spike conversion ...');
    numTst = length(PtnCellTst);
    for i = 1:numTst
        for j = 1:size(PtnCellTst{i}.AllVec,2)
            vec = PtnCellTst{i}.AllVec(:,j);
            vec(vec>maxResp) = maxResp;
            PtnCellTst{i}.AllVec(:,j) = 1- vec / maxResp ;
        end
    end
    save (file, 'PtnCellTst')
end
load (file)




%%

file = [folder0,'/','TrainedWt.mat'];
load (file)


file = [folder1,'/','RawResults.mat'];
if ( (reRunAll) || (reSplitDataset) || ~exist(file,'file') )
    timedLog('classification AllTstIn1 ...');
    IsTraining = 0;
%     [TrainedWt, RR_tst, RR_tst_TP, RR_tst_TN]=edTempotronClassify(TrainedWt, IsTraining, 'testing set',folder1);
    [TrainedWt, RR_tst, RR_tst_TP, RR_tst_TN]=edTempotronClassify_callMEX(TrainedWt, IsTraining, 'testing set',folder1); 
    save (file, 'RR_tst', 'RR_tst_TP', 'RR_tst_TN');
end  
load (file)
fprintf('======================================================\n');
fprintf('Raw correct Rate (i.e. rate of correctly responded neurons for all slices):\n');
fprintf('\t AllTstIn1 : total succ rate %.2f %%, \t true positive %.2f %%, \t true negative %.2f %% \n',  ...
    RR_tst*100, RR_tst_TP*100, RR_tst_TN*100);



% analyze results (majority voting) to make final decisions
file = [folder1,'/','FinalResults.mat'];
if ( (reRunAll) || (reSplitDataset) || ~exist(file,'file') )
    [ FinalCorrRate_Tst, decision_Tst, freqsCell_Tst, ConfMatrix_Tst ] = ...
        AnalyzeResults_MajorityVoting_AllTstIn1( nOutputs, nNeuronPerOutput, folder1 );
    save (file, 'FinalCorrRate_Tst', 'decision_Tst', 'freqsCell_Tst', 'ConfMatrix_Tst')
end
load (file)
fprintf('======================================================\n');
fprintf('Final correct Rate (i.e. succ rate of categorization, by majority voting):\n');
fprintf('\t AllTstIn1 :  %.2f %% \n',  FinalCorrRate_Tst*100);

fprintf('Confusion Matrix : column index as output, row index as label\n');
disp(ConfMatrix_Tst)


% draw the figure, ground truth and algorithm decisions
file = [folder1, '/', 'PtnCellTst_out.mat'];
load (file, 'PtnCellTst')
tdecision = [PtnCellTst{1}.Time_Chnl_Lbl(1,:)', decision_Tst{1}'];
h=figure; 
plot(tlabel(:,1), tlabel(:,2),'-b', tdecision(:,1),tdecision(:,2),'.r');
legend('ground truth', 'decision made');
xlabel('time(ns)'); ylabel('class'); 
title(['AllTstIn1CIN, corr rate = ',num2str(FinalCorrRate_Tst*100,'%.2f'),' %']);


