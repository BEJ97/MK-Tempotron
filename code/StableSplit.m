function [PtnCellTrn, PtnCellTst, indTrn, indTst, indNaN] = StableSplit(PtnCell,ratio_Trn,nGroup)
% RANDSPLIT Split the dataset into two parts (train and test) randomly,
% the ratio of training samples to the whole set is specified by ratio_Trn.

timedLog('randomly splitting the dataset...');

N = length(PtnCell);
%N = 1000;
Labels = zeros(1,N);   %改 原Labels = zeros(1,N); 
numTotSlices = 0;
for i = 1:N        %改原for i = 1:N
    B = PtnCell{i}.Time_Chnl_Lbl;
    numTotSlices = numTotSlices + size(B,2);
    if ~isempty(B),
        Labels(i) = B(3,1);
    else
        Labels(i) = NaN;
    end
end

indTrn = zeros(1, round(N*ratio_Trn));
indTst = zeros(1, N-round(N*ratio_Trn));
count1 = 0;
count2 = 0;
for i = 1:nGroup
    ind = find(Labels==(i-1));
    n0 = length(ind);
    
	%%%
	ind2 = zeros(1,n0);
	for j = 1 : n0
		ind2(1, j) = j;
	end
	%%%%
    n1 = round(n0*ratio_Trn);
    n2 = n0-n1;
    indTrn(count1+1:count1+n1) = ind(ind2(1:n1));
    indTst(count2+1:count2+n2) = ind(ind2(n1+1:n0));
    count1 = count1+n1;
    count2 = count2+n2;
end

	% i = 1
    % ind = find(Labels==(i-1));
    % n0 = length(ind);
    % ind2 = randperm(n0);
    % n1 = round(n0*ratio_Trn);
    % n2 = n0-n1;
    % indTrn(count1+1:count1+n1) = ind(ind2(1:n1));
    % indTst(count2+1:count2+n2) = ind(ind2(n1+1:n0));
    % count1 = count1+n1;
    % count2 = count2+n2;


indTrn = indTrn(1:count1);
indTst = indTst(1:count2);
indNaN = find(isnan(Labels));

numTrn = length(indTrn);
numTst = length(indTst);
numNaN = length(indNaN);

indTrn = indTrn(randperm(numTrn));
indTst = indTst(randperm(numTst));



PtnCellTrn = PtnCell(indTrn);
PtnCellTst = PtnCell(indTst);


    
    
    


