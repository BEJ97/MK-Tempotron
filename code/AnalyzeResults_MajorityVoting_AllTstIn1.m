function [ CorrRate, decision, freqsCell, ConfMatrix ] = ...
    AnalyzeResults_MajorityVoting_AllTstIn1( nOutput, nNeuronPerOutput, datafolder )

if ~exist('datafolder','var')
    datafolder = pwd;
end

load ([datafolder,'/','PtnCellTst_out.mat'])



%-------------- correct rate ---------------

N = length(PtnCellTst);
% decision = zeros(1,N);
decision = cell(1,N);
freqsCell = cell(1,N);
numCorrect = 0;
numTotSlice = 0;

for i = 1:N

    Out = PtnCellTst{i}.Out;
    nSlices = size(Out,2);
    freqs = zeros(nOutput,nSlices);
    
    
    for indSlice = 1:nSlices
        for neuron = 1:nOutput
            for indNeuronPerOutput = 1:nNeuronPerOutput

                indTotNeuronOutput = ((neuron-1)*nNeuronPerOutput+indNeuronPerOutput);
                out = Out(indTotNeuronOutput,indSlice);

                if out==1
                    freqs(neuron,indSlice) = freqs(neuron,indSlice) +1;
                end
            end
        end
    end
    
    freqsCell{i} = freqs;

    decision{i} = zeros(1,nSlices);
    for indSlice = 1:nSlices
        freqs1 = freqs(:,indSlice);
        decision{i}(indSlice) = find(freqs1==max(freqs1),1) - 1;
        lbl    = PtnCellTst{i}.Time_Chnl_Lbl(3,indSlice);
        numCorrect = numCorrect + double(decision{i}(indSlice)==lbl(1));
        numTotSlice = numTotSlice + 1;
    end

end

CorrRate = numCorrect/numTotSlice;




%-------------- conf matrix --------------

nGrp = nOutput;
ConfMatrix = zeros(nGrp+1,nGrp+1);  % target class (actual class), output class (predicted class)

N = length(PtnCellTst);
for i = 1:N
    Out = PtnCellTst{i}.Out;
    nSlices = size(Out,2);
    for indSlice = 1:nSlices
        y = decision{i}(indSlice);                      % output class (predicted class)
        t = PtnCellTst{i}.Time_Chnl_Lbl(3,indSlice);    % target class (actual class)
        ConfMatrix(t+1,y+1) = ConfMatrix(t+1,y+1) +1;
    end
end
for i = 1:nGrp
    ConfMatrix(i,nGrp+1) = ConfMatrix(i,i)/sum(ConfMatrix(i,1:nGrp));
end
for j = 1:nGrp
    ConfMatrix(nGrp+1,j) = ConfMatrix(j,j)/sum(ConfMatrix(1:nGrp,j));
end
ConfMatrix(nGrp+1,nGrp+1) = sum(diag(ConfMatrix(1:nGrp,1:nGrp))) / sum(sum(ConfMatrix(1:nGrp,1:nGrp)));


end




