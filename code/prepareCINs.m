function [NumCINs] = prepareCINs(datafolder)

if ~exist('datafolder','var')
    datafolder = pwd;
end


load ('posture_2timesIn1CIN_32x32.mat')%改

for i = 1:length(CINs)
    CINs{i}(:,1) = CINs{i}(:,1)*1e3;        % us --> ns
end

NumCINs = length(CINs);

save ([datafolder,'/','posture.mat'], 'CINs', 'Labels','NumCINs')

end