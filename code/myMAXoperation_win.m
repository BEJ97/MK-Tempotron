function [ stateCell ] = myMAXoperation_win( stateCell,filtersize,iEvt, event_in,sqfilters,post_proc_th)
% myMAXoperation_win s1->c1: max over local
% do max for a window region only.

[numScales, numRot] = size(stateCell);

x = event_in(1,4);
y = event_in(1,5);

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
        
        c = floor(size(sqfilters{count},1)/2);            
        winx = max(1,x-c):min(M,x+c); 
        winy = max(1,y-c):min(N,y+c); 

        % only need to do max for pixels in this window
        for rr= winx           
            for cc= winy   
                r1=max(1,rr-sz); r2=min(M,rr+sz);
                c1=max(1,cc-sz); c2=min(N,cc+sz);
                if tmp(rr,cc)<post_proc_th 
                    tmp2(rr,cc) = 0;    % thresholding to remove small values.
                elseif tmp(rr,cc)~=max(max( tmp(r1:r2,c1:c2) ))
                    tmp2(rr,cc) = 0;    % kill non-max.
                else % MAX surviving

                end
            end
        end
        
        stateCell{idxScale,idxRot}.J_max(winx,winy)=tmp2(winx,winy);
    
        stateCell{idxScale,idxRot}.time_maxCheck_b(winx,winy) = event_in(1,1)*ones(length(winx),length(winy));
    end
end




