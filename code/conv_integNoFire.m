function [state]=conv_integNoFire(event_in,params,state,CONV_refresh_window_only)
% convolution, leaky integration. no fire.

t   = event_in(1,1);
x   = event_in(1,4);                % row addr
y   = event_in(1,5);                % col addr
sgn = event_in(1,6);                % sign

size1      = state.size1;
size2      = state.size2;
J          = state.J;               % the response map after convolution
time_in_b  = state.time_in_b;       % a matrix that stores the time of last input event for each neuron. -inf by default

cteloss     = params.cteloss;       % leakage rate.
s           = params.s;             % filter, convolution mask

z1 = params.zs(1);                  % These four parameters specify the window of operation.
z2 = params.zs(2);                  % e.g. z1=-1, z2=1, z3=-1,  z4=1.
z3 = params.zs(3);                  % The window is ( x+z1: x+z2, y+z3: y+z4 )
z4 = params.zs(4);                  % The window is also limited by the boundaries.               
winx = max(1,x+z1):min(size1,x+z2); % so, the final window (winx, winy) is 
winy = max(1,y+z3):min(size2,y+z4); % (max(1,x+z1):min(size1,x+z2), max(1,y+z3):min(size2,y+z4))

% a larger region than (winx,winy). 
% Calc leakage for region (winx2,winy2). 
% But adding the kernle only to region (winx,winy).
winx2 = max(1,x+2*z1):min(size1,x+2*z2);
winy2 = max(1,y+2*z3):min(size2,y+2*z4);


% effective mask, to deal with boundary problems -------------------
ind11 = max(1,2-x-z1);
ind21 = max(1,2-y-z3);
ind12 = min(size(s,1),(size(J,1)+size(s,1)-(x+z2)));
ind22 = min(size(s,2),(size(J,2)+size(s,2)-(y+z4)));
s2    = s(ind11:ind12,ind21:ind22); 


% forgetting mechanism
if CONV_refresh_window_only,
    if cteloss~=0
        leakage  = cteloss* ( t-time_in_b(winx2,winy2) );
        J(winx2,winy2) = sign(J(winx2,winy2)).* max( abs(J(winx2,winy2))- leakage, 0);
    end   
    % update the time_in_b.
    state.time_in_b(winx2,winy2) = t*ones(length(winx2),length(winy2)); 
else
    if cteloss~=0
        leakage  = cteloss*( t-time_in_b );
        J = sign(J).*( abs(J) - min( abs(J), abs(leakage) ) ); 
    end   
    % update the time_in_b.
    state.time_in_b = t*ones(size1,size2); 
end

                                    
if ~isempty(s2)
    % add kernel to the feature map.
    J(winx, winy) = J(winx, winy)+ sgn*s2;     
                                             
    state.J = J;    
    
end

