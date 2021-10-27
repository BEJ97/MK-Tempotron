function [ SpkTimings,Addresses,TimeChnlLbl,NeuronsTgt, NeuronsOut, ...
    idx_convert, hasFieldTgt, hasFieldOut] = ...
    ConvertPtnCell2Arrays( PtnCell,nOutputs,nNeuronPerOutput )

% ---------- PtnCell to Arrays:
% PtnCell    training   testing     
% SpkTimings 100x60000, 100x10000   double
% Addresses  100x60000, 100x10000   double
% Labels     1x60000,   1x10000     double
% NeuronsTgt 100x60000  100x10000   logical
% NeuronsOut 100x60000  100x10000   logical

nImages = length(PtnCell);
nSurvC1 = size(PtnCell{1}.AllVec,1);

SpkTimings = zeros(nSurvC1, nImages);
Addresses  = zeros(nSurvC1, nImages);
TimeChnlLbl= zeros(3, nImages);
NeuronsTgt = false(nOutputs*nNeuronPerOutput, nImages);
NeuronsOut = false(nOutputs*nNeuronPerOutput, nImages);

idx_convert = zeros(2, nImages); % record iImage and iSlice, for reverse conversion.

cnt = 0;
for iImage = 1:nImages
	nSlices = size(PtnCell{iImage}.AllVec,2);
	cols = (cnt+1) : (cnt+nSlices) ;
	
	SpkTimings(:, cols) = PtnCell{iImage}.AllVec;
	Addresses(:, cols)  = PtnCell{iImage}.AllAddr;
	TimeChnlLbl(:,cols) = PtnCell{iImage}.Time_Chnl_Lbl(1:3,:);
    
    hasFieldTgt = isfield(PtnCell{iImage},'Tgt');
    if hasFieldTgt
        NeuronsTgt(:, cols) = PtnCell{iImage}.Tgt;
    end
    hasFieldOut = isfield(PtnCell{iImage},'Out');
    if hasFieldOut
        NeuronsOut(:, cols) = PtnCell{iImage}.Out;
    end
    
    tmp = [];
    for iSlice = 1:nSlices
        tmp2 = [iImage; iSlice];
        tmp = [tmp, tmp2];
    end
    idx_convert(:,cols) = tmp;
	
	cnt = cnt+nSlices;
end

% The final column size may be larger than nImages!
% Correct the size to avoid possible errors in later MEX function!
if ~hasFieldTgt
    NeuronsTgt = false(nOutputs*nNeuronPerOutput, size(SpkTimings,2));
end
if ~hasFieldOut
    NeuronsOut = false(nOutputs*nNeuronPerOutput, size(SpkTimings,2));
end


% ----------------


end

