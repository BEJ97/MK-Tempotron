function [ stateCell ] = myMAXoperation( stateCell,filtersize,post_proc_th)
% myMAXoperation s1->c1: max over local
% do max for whole map!

[numScales, numRot] = size(stateCell);

% ---------------------------------------------------
% "max" over local neighborhood
count = 0;
for idxScale = 1:numScales
    for idxRot = 1:numRot
        count = count +1;
        tmp= abs( stateCell{idxScale,idxRot}.J );  % S1.
        tmp2 = tmp;
        sz = floor(filtersize(count)/2);
        [M,N] = size(tmp);

        % do max for whole map!
        for rr=1:M
            for cc=1:N
                r1=max(1,rr-sz); r2=min(M,rr+sz);
                c1=max(1,cc-sz); c2=min(N,cc+sz);         
                if tmp(rr,cc)<post_proc_th 
                    tmp2(rr,cc) = 0;    % thresholding to remove small values.
                elseif tmp(rr,cc)~=max(max( tmp(r1:r2,c1:c2) ))
                    tmp2(rr,cc) = 0;    % kill non-max.
                end
            end
        end

        stateCell{idxScale,idxRot}.J_max=tmp2;
        
    end
end




