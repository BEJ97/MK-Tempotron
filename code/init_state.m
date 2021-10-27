function [ state ] = init_state(size1,size2  )
%INIT_STATE set initial values of parameters in "state" struct.

if ~exist('size1','var')
    size1 = 32;
end
if ~exist('size2','var')
    size2 = 32;
end

state.size1 = size1;
state.size2 = size2;
state.J = zeros(size1,size2);           % S1 map
state.J_max = zeros(size1,size2);       % C1 map
state.time_in_b  = -inf*ones(size1,size2);   
state.time_maxCheck_b  = -inf*ones(size1,size2); 

end

