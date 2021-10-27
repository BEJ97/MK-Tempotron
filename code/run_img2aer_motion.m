function run_img2aer_motion(d, ratio, datafolder,nRepeat)

if ~exist('datafolder','var')
    datafolder = pwd;
end

if ~exist('nRepeat','var'),
    nRepeat = 10;
end

timedLog('running img2aer...');

load MNIST_blackBg_60000train_10000test

Images = [TrainImages,TestImages];
Labels = [TrainLabels,TestLabels];


n = length(Labels);
CINs = cell(1,n);
for i = 1:n
    I = reshape(Images(:,i),28,28)';
    th = 255*0.4;
    I = double(I>th)*255;
    Images(:,i) = reshape(I',28*28,1);
    CIN = img2aer_motion(I,d,ratio,nRepeat);
    CINs{i} = CIN;
end

NumCINs = n;

save ([datafolder,'/','AER_MNISTbRand_500_alt'], 'Images', 'CINs', 'Labels','NumCINs')




