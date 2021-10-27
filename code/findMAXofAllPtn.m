function [maxi] = findMAXofAllPtn(PtnCell)

N = length(PtnCell);
maxi = -inf;

for i = 1:N
    AllVec = PtnCell{i}.AllVec;
    if isempty(AllVec)
        disp(i);
    end
    maxi = max(maxi,max(AllVec(:)));
end


end
