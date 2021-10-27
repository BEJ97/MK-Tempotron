function [ CorrRate_Trn, CorrRate_Tst, decision_Trn, decision_Tst, freqsCell_Trn, freqsCell_Tst ] = ...
    AnalyzeResults_MajorityVoting( nOutput, nNeuronPerOutput, datafolder )

if ~exist('datafolder','var')
    datafolder = pwd;
end

load ([datafolder,'/','PtnCellTrn_out.mat'])
[CorrRate_Trn, decision_Trn, freqsCell_Trn] = AnalyzeResults_MajorityVoting_Core (PtnCellTrn, ...
    nOutput,nNeuronPerOutput);

load ([datafolder,'/','PtnCellTst_out.mat'])
[CorrRate_Tst, decision_Tst, freqsCell_Tst] = AnalyzeResults_MajorityVoting_Core (PtnCellTst, ...
    nOutput,nNeuronPerOutput);

end




function [CorrRate, decision,freqsCell] = AnalyzeResults_MajorityVoting_Core (PtnCellTst, ...
    nOutput,nNeuronPerOutput)

N = length(PtnCellTst);
% decision = zeros(1,N);
decision = cell(1,N);
freqsCell = cell(1,N);
numCorrect = 0;
numTotSlice = 0;

for i = 1:N    
    
    lbl    = PtnCellTst{i}.Time_Chnl_Lbl(3,1);

    Out = PtnCellTst{i}.Out;
    nSlices = size(Out,2);
    freqs = zeros(nOutput,nSlices);
    
    
    for indSlice = 1:nSlices
        for neuron = 1:nOutput   %��ԭû%
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
        numCorrect = numCorrect + double(decision{i}(indSlice)==lbl(1));
        numTotSlice = numTotSlice + 1;
    end


end

CorrRate = numCorrect/numTotSlice;

end

