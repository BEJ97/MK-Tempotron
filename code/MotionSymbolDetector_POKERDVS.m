function [ Timings, V, t ] = MotionSymbolDetector_POKERDVS( CIN, ...
    V0,dt,use_single_exponential,lut1,lut2, THs, label, nGroup, size1, plotfig)
%MotionSymbolDetector
% Inputs: 
%     CIN, address events
%     V0, dt, use_single_exponential, lut1,lut2, : related to exp decay lookup table 
%     THs: several different thresholds, use min(THs) to remove small peaks caused by noise
%     label: the label for current CIN.
%     nGroup: number of classes.
%     size1, input spatial resolution. eg. 32
%     plotfig: if 1, plot the integration figure.
% Output: 
%     Timings, 1st row: time; 2nd row: channel (0~9); 3rd row: label;
%     V: potential of the leaky integration neuron
%     t: time stamps of input address events.  

if ~exist('plotfig','var'),
    plotfig = 0;
end

t = CIN(:,1);
[ V ]= edInteg (t,V0,dt,use_single_exponential,lut1,lut2) ;

if (plotfig),
    figure,plot(t,V,'-k'); xlabel('t(ns)'); ylabel('integ potential'); xlim([t(1),t(end)]);
    drawnow;
end

Timings = zeros(3, nGroup*10);  
% 1st row: time; 2nd row: channel (0~9); 3rd row: label; 
count = 0;

th = THs(label+1);           %改       th = min(THs); 2020/7/20

timeRadius = 700;  % 15ms  改 原4e6

isi = mean( t(2:end)-t(1:end-1) );
Radius = 700;

[ tpeak ] = FindPeak(V,t,Radius,th);

[ tExact ] = getExactEvtTime(t, tpeak);

tmp = [tExact; 1*ones(size(tExact)); label*ones(size(tExact))]; 
if ~isempty(tExact),
    Timings(:,count+1:count+length(tExact)) = tmp;
end
count = count+ length(tExact);

Timings = Timings(:,1:count);

[~, IX] = sort(Timings(1,:));
Timings = Timings(:,IX);



