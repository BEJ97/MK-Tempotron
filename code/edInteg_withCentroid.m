function [V, centroid,ts] = edInteg_withCentroid (CIN,size1,...
    V0,dt,use_single_exponential,lut1,lut2)
% leaky integration neuron (in Motion Symbol Detector) as well as centroid calculation.

ts = CIN(:,1);
NumEvt = length(ts);

t_last = -1;
K1 = zeros(size1,size1);
K2 = zeros(size1,size1);

V  = zeros(NumEvt,1);       % total potential
xc = zeros(NumEvt,1);       % centroid x
yc = zeros(NumEvt,1);       % centroid y

for i = 1:NumEvt
    ti = ts(i);
    xi = CIN(i,4);
    yi = CIN(i,5);
    
    delta_t = ti-t_last;
    
    lut_addr = round(delta_t /dt)+1;
    if lut_addr<=length(lut1)
        Sc1 = lut1(lut_addr);
    else
        Sc1 = 0;
    end
    K1 = Sc1*K1;
    K1(xi,yi) = K1(xi,yi) + V0;
    
    if ~use_single_exponential,
        if lut_addr<=length(lut2)
            Sc2 = lut2(lut_addr);
        else
            Sc2 = 0;
        end
        K2 = Sc2*K2;
        K2(xi,yi) = K2(xi,yi) + V0;
        K = K1-K2;
    else
        K = K1; 
    end
 
    V(i) = sum(sum(K));             % total potential. (i.e. the potential in Motion symbol detector)
    
    denom = V(i)+ ~V(i);            % change 0 with 1.
    xc(i) = ( sum( (1:size1) *K ) / denom );
    yc(i) = ( sum( (1:size1) *K') / denom );
 
    t_last = ti;
end

xc = ceil( xc + ~xc *size1/2 );
yc = ceil( yc + ~yc *size1/2 );

centroid = [xc';yc'];

