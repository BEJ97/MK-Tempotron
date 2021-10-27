function [ paramCell,stateCell,sqfilters,filtersize,numScales,numRot,numFilters] = ...
    init_paramCell_stateCell( minFS,maxFS,sSFS,rot,circularRField,useSmooth, cteloss,size1,size2 )

[sqfilters,filtersize,filters]= zb_filters2(minFS,maxFS,sSFS,rot,circularRField,0,useSmooth);

RF_siz = minFS:sSFS:maxFS;
numScales = length(RF_siz);
numRot = length(rot);
numFilters = numScales*numRot;

params.cteloss = cteloss;            % lose how much per time unit (nanosecond) 

paramCell = cell(numScales,numRot);
for i = 1:numScales
    for j = 1:numRot
        idx = (i-1)*numRot+j ;
        params.s = sqfilters{idx};          % filter,2d
        c = floor(size(sqfilters{idx},1)/2);
        params.zs = [-c c -c c];            % operation window
        paramCell{i,j} = params;
    end
end


state = init_state(size1,size2);

stateCell = cell(numScales,numRot);
for i = 1:numScales
    for j = 1:numRot
        stateCell{i,j} = state;
    end
end


end

