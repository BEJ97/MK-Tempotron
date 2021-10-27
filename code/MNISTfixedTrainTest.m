function [PtnCellTrn, PtnCellTst, indTrn, indTst,indNaN] = MNISTfixedTrainTest(PtnCell,numTrainingImgs)
% The first 60000 for training, the other 10000 for testing.

ntest = 70000-numTrainingImgs;
timedLog(sprintf('MNIST, %d for training, %d for testing...',numTrainingImgs,ntest));

indTrn = 1:numTrainingImgs;
indTst = numTrainingImgs+1:70000;

PtnCellTrn = PtnCell(indTrn);
PtnCellTst = PtnCell(indTst);

indNaN = [];


    
    
    


