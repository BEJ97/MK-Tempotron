function [ PtnCell ] = ConvertArrays2PtnCell(SpkTimings,Addresses,...
    TimeChnlLbl,NeuronsTgt,NeuronsOut,idx_convert,hasFieldTgt,hasFieldOut)

numTotSlices = size(SpkTimings,2);
nImages = max(idx_convert(1,:));
PtnCell = cell(1,nImages);
% ---------- Arrays to PtnCell :
for pp = 1: numTotSlices
    iImage = idx_convert(1,pp);
    iSlice = idx_convert(2,pp);
    PtnCell{iImage}.AllVec(:,iSlice) = SpkTimings(:,pp);
    PtnCell{iImage}.AllAddr(:,iSlice)= Addresses(:,pp);
    PtnCell{iImage}.Time_Chnl_Lbl(:,iSlice) = TimeChnlLbl(:,pp);
    if hasFieldTgt
        PtnCell{iImage}.Tgt(:,iSlice) = NeuronsTgt(:, pp);
    end
    if hasFieldOut
        PtnCell{iImage}.Out(:,iSlice) = NeuronsOut(:, pp);
    end
end
% ----------------


end

