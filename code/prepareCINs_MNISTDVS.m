function [NumCINs] = prepareCINs_MNISTDVS(datafolder,sourceMatFile)

if ~exist('datafolder','var')
    datafolder = pwd;
end
if ~exist('sourceMatFile','var')
    sourceMatFile = 'MNIST_DVS_28x28_10000_full.mat';
end


load (sourceMatFile)

for i = 1:length(CINs)
    CINs{i}(:,1) = CINs{i}(:,1)*1e3;        % us --> ns
end

NumCINs = length(CINs);

save ([datafolder,'/','MNISTDVS.mat'], 'CINs', 'Labels','NumCINs')

end