function [ Vec ] = myRasterScan(C1)
% MYRASTERSCAN reshape C1 cell to a column vector.

[numScales, numRot] = size(C1);
Vec = [];

for i = 1:numScales
    for j = 1:numRot
        [M,N] = size(C1{i,j});
        tmp = reshape(C1{i,j}', M*N,1);
        Vec = [Vec; tmp];
    end
end


end
