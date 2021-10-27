function [NumCINs] = prepareCINs_CIFAR10_DVS(datafolder,sourceMatFile)

if ~exist('datafolder','var')
    datafolder = pwd;
end
if ~exist('sourceMatFile','var')
    sourceMatFile = 'CIFAR10-DVS0.mat';
end


load (sourceMatFile)

for i = 1:length(CINs)
    CINs{i}(:,1) = CINs{i}(:,1)*1e6;        % ms --> ns
end

NumCINs = length(CINs);

save ([datafolder,'/','CIFAR10-DVS.mat'], 'CINs', 'Labels','NumCINs')

end