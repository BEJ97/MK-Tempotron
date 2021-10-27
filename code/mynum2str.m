function [ str ] = mynum2str( x,f )
%MYNUM2STR 
% MATLAB windows version and unix version gives different results for num2str(5e-5).
% '5e-005' (for windows version), '5e-05' for unix version.
% This function is to modify unix version's result to make it consistent
% with that of windows version.

if nargin<1
    str = num2str();
elseif nargin<2
    str = num2str(x);
    if length(str)>3
        if str(end-3)=='e'
            str = [str(1:end-2) '0' str(end-1:end)];
        end
    end
elseif ( ischar(f) && (f(end)=='e') )
    str = num2str(x,f);
    if str(end-2)=='-' || str(end-2)=='+'
        str = [str(1:end-2) '0' str(end-1:end)];
    end
else
    str = num2str(x,f);
end


end

