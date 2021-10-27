function [ Timings, V, t ] = MotionSymbolDetector_withCentroid_AllTstIn1( CIN, ...
    V0,dt,use_single_exponential,lut1,lut2, THs, tlabel, nGroup, size1, plotfig)
%MotionSymbolDetector_withCentroid_AllTstIn1
% Inputs: 
%     CIN, address events
%     V0, dt, use_single_exponential, lut1,lut2, : related to exp decay lookup table 
%     THs: several different thresholds, use min(THs) to remove small peaks caused by noise
%     tlabel: ground truth of the event stream. (time, label)
%     nGroup: number of classes.
%     size1, input spatial resolution. eg. 32
%     plotfig: if 1, plot the integration figure.
% Output: 
%     Timings, 1st row: time; 2nd row: channel (0~9); 3rd row: label;
%               4th row: centroid row addr; 5th row: centroid col addr
%     V: potential of the leaky integration neuron
%     t: time stamps of input address events.

if ~exist('plotfig','var'),
    plotfig = 0;
end

[V, centroid,t] = edInteg_withCentroid (CIN,size1,...
    V0,dt,use_single_exponential,lut1,lut2) ;
% centroid: 2xN matrix. 1st row: centroid row addr, 2nd row: centroid column addr.

if (plotfig),
    figure,
    subplot(211),plot(t,V,'-k'); xlabel('t(ns)'); ylabel('integ potential'); xlim([t(1),t(end)]);
    subplot(212),plot(t,centroid(1,:),'.r',t,centroid(2,:),'.b'); legend('row','col'); xlabel('t(ns)'); ylabel('centroid addr');  xlim([t(1),t(end)]);
    drawnow;
end


Timings = zeros(5, nGroup*10);  
% 1st row: time; 2nd row: channel (0~9); 3rd row: label; 4th row: centroid row addr; 5th row: centroid col addr
count = 0;

th = min(THs);

timeRadius = 1e9; %unit: ns, 1e9ns=1s
% timeRadius = 5e8;   % 500ms
% timeRadius = 2e8; % 200ms

isi = mean( t(2:end)-t(1:end-1) );
Radius = ceil(timeRadius/isi);

[ tpeak,centroid_atPeak ]=FindPeak_withCentroid(V,t,centroid,Radius,th) ;

[ tExact ] = getExactEvtTime(t, tpeak);

labels = zeros(1,length(tExact));
for i = 1:length(tExact)
    ts = tlabel(:,1);
    e = find( ( ts(1:end-1)<=tExact(i) ) & ( ts(2:end)>=tExact(i) ) );
    labels(i) = tlabel(e(1),2);
end

tmp = [tExact; 1*ones(size(tExact)); labels; centroid_atPeak]; 
if ~isempty(tExact),
    Timings(:,count+1:count+length(tExact)) = tmp;
end
count = count+ length(tExact);

Timings = Timings(:,1:count);

[~, IX] = sort(Timings(1,:));
Timings = Timings(:,IX);


