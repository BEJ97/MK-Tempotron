function [CIN_sorted] = img2aer_motion(I,d,ratio,nRepeat)
% binary image --> motion aer events.
% inputs:
%   I:  binary image
%   d:  salt & pepper noise density
%   ratio:  the ratio of pixels that generate events to all white pixels
% outputs:
%   CIN_sorted: Nx6 matrix
%       [nanosecond timestamp,  0, -1, x address, y address,  polarity]


% clear;clc;
% load MNIST_blackBg_whiteFg
% 
% I = double(TestImages{1});
% I = double(I>0)*255;
% imshow(I);

if ~exist('nRepeat','var'),
    nRepeat = 10;
end


% tic

I = double(I);
[M,N]=size(I);


% T = 2e-3;             % unit: second. T=2ms
%n = 10;                 % each foreground pixel generates about 10 events.
n = nRepeat;		
ISI = 1e-7;             % 100ns between two adjacent events
% step = 1/4*M*N *ISI* 10;% interval of repeatition (same address)
step = 1/4*M*N *ISI* 5; % interval of repeatition (same address)

% d = 0.05;               % salt & pepper noise density
% ratio = 0.9 ;           % the ratio of active pixels that generate events.


CIN = zeros(M*N*n,3);  
count = 0;
for iRepeat = 1:n
    J = imnoise(I,'salt & pepper',d);
    [row,col] = find(J~=0);
    
    %ind = randperm(length(row), round(ratio*length(row)) );
    ind = randperm(length(row));
    ind = ind( 1:round(ratio*length(row)) ) ;
    row = row(ind);
    col = col(ind);
    
    ts = (iRepeat-1)*step + (0:length(row)-1)'*ISI;
    temp = [ts, row, col];
       
    numEvt = length(ind);
    CIN(count+1:count+numEvt,:) = temp;
    count = count + numEvt;
end
    
    
CIN = CIN(1:count,:);
% [B,IX]=sort(CIN(:,1));
% CIN_sorted = CIN(IX,:);
CIN_sorted = CIN;

CIN_sorted = [round(CIN_sorted(:,1)*1e9), zeros(count,1),-ones(count,1), ...
    CIN_sorted(:,2), CIN_sorted(:,3), ones(count,1)];
% nano second
% 
% toc


% %%
% A = CIN_sorted;
% figure,scatter3(A(:,4),A(:,5),A(:,1));xlim([1 28]);ylim([1 28]);
% figure,scatter3(A(:,5),A(:,4),A(:,1));xlim([1 28]);ylim([1 28]);axis ij;
% figure,scatter3(A(:,5),A(:,4),A(:,1)/1e3,'.');xlim([1 28]);ylim([1 28]);axis ij;


            