function [NumCINs] = prepareCINs_N_CARS(datafolder,sourceMatFile)

if ~exist('datafolder','var')
    datafolder = pwd;
end
if ~exist('sourceMatFile','var')
    sourceMatFile = 'N-CARS_all.mat';
end


load (sourceMatFile)

for i = 1:length(CINs)
    CINs{i}(:,1) = CINs{i}(:,1)*1e6;        % ms --> ns
end

NumCINs = length(CINs);

save ([datafolder,'/','N-CARS.mat'], 'CINs', 'Labels','NumCINs')

end