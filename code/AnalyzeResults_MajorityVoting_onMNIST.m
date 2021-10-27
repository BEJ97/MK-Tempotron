function [ CorrRate_Trn, CorrRate_Tst, decision_Trn, decision_Tst, freqsCell_Trn, freqsCell_Tst ] = ...
    AnalyzeResults_MajorityVoting_onMNIST( nOutput, nNeuronPerOutput, datafolder )
%AnalyzeResults_MajorityVoting_onMNIST

if ~exist('datafolder','var')
    datafolder = pwd;
end

load ([datafolder,'/','PtnCellTrn_out.mat'])
[CorrRate_Trn, decision_Trn, freqsCell_Trn] = AnalyzeResults_MajorityVoting_onMNIST_Core (PtnCellTrn, ...
    nOutput,nNeuronPerOutput);

load ([datafolder,'/','PtnCellTst_out.mat'])
[CorrRate_Tst, decision_Tst, freqsCell_Tst] = AnalyzeResults_MajorityVoting_onMNIST_Core (PtnCellTst, ...
    nOutput,nNeuronPerOutput);

end




function [CorrRate, decision,freqsCell] = AnalyzeResults_MajorityVoting_onMNIST_Core (PtnCellTst,...
    nOutput,nNeuronPerOutput)

N = length(PtnCellTst);
decision = zeros(1,N);
freqsCell = cell(1,N);
numCorrect = 0;

for i = 1:N
    
    lbl    = PtnCellTst{i}.Time_Chnl_Lbl(3,1);

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
    
    freqs1 = sum(freqs,2);
    decision(i) = find(freqs1==max(freqs1),1) - 1;
    numCorrect = numCorrect + double(decision(i)==lbl(1));
end

CorrRate = numCorrect/N;

end

