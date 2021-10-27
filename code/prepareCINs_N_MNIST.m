function [NumCINs] = prepareCINs_N_MNIST(datafolder,sourceMatFile)

if ~exist('datafolder','var')
    datafolder = pwd;
end
if ~exist('sourceMatFile','var')
    sourceMatFile = 'N-MNIST_all.mat';
end


load (sourceMatFile)

for i = 1:length(CINs)
    CINs{i}(:,1) = CINs{i}(:,1)*1e3;        % us --> ns
end

NumCINs = length(CINs);

save ([datafolder,'/','NMNIST.mat'], 'CINs', 'Labels','NumCINs')

end