function [ ] = myC1respToSpk (maxResp, datafolder)
% myC1respToSpk

if ~exist('datafolder','var')
    datafolder = pwd;
end

timedLog('c1 resp --> TFS spikes...');

load ([datafolder,'/','PtnCellTrnTst_raw'], 'PtnCellTrn', 'PtnCellTst', 'indTrn', 'indTst', 'indNaN')

numTrn = length(PtnCellTrn);
numTst = length(PtnCellTst);

for i = 1:numTrn
    for j = 1:size(PtnCellTrn{i}.AllVec,2)
        vec = PtnCellTrn{i}.AllVec(:,j);
        vec(vec>maxResp) = maxResp;
        PtnCellTrn{i}.AllVec(:,j) = 1- vec/maxResp ;%+ rand(1,1)/40 ;  %改
    end
end

for i = 1:numTst
    for j = 1:size(PtnCellTst{i}.AllVec,2)
        vec = PtnCellTst{i}.AllVec(:,j);
        vec(vec>maxResp) = maxResp;
        PtnCellTst{i}.AllVec(:,j) = 1- vec/maxResp ;%+ rand(1,1)/40 ;   %改
    end
end

save ([datafolder,'/','PtnCellTrnTst_spk'], 'PtnCellTrn', 'PtnCellTst', 'indTrn', 'indTst', 'indNaN')

