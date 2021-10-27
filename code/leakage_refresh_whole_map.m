function [state]= leakage_refresh_whole_map(params,state,timeact)

size1      = state.size1;
size2      = state.size2;
J          = state.J;               % the response map after convolution. S1
J_max      = state.J_max;           % C1.
time_in_b  = state.time_in_b;       % a matrix that stores the time of last input event for each neuron. -inf by default
time_maxCheck_b = state.time_maxCheck_b;


cteloss     = params.cteloss;       % leakage rate.

% forgetting mechanism, whole map refreshing
if cteloss~=0
	leakage = cteloss*(timeact-time_in_b);
    J = sign(J).*( abs(J) - min( abs(J), abs(leakage) ) ); 

    leakage2 = cteloss*(timeact-time_maxCheck_b);
    J_max = sign(J_max).*( abs(J_max) - min( abs(J_max), abs(leakage2) ) ); 
end   

% update the time_in_b.
state.time_in_b = timeact * ones(size1,size2); 
state.time_maxCheck_b = timeact * ones(size1,size2); 

state.J = J;
state.J_max = J_max;


