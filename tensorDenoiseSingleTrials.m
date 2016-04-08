%% Script: single trial shenanigans 

%%
load('/Users/Jeff/Documents/MATLAB/TensorDenoising/mat-files/DataLaraB.mat')
DataT = Data;
DataT.Y = DataT.Y(:,:,:);
DataT.Ys = DataT.Ys(:,:,:);
DataT.Ysm = DataT.Ysm(:,:,:);
[n,t,c,r] = size(Data.Y);

%% get 'true' rank
[~, mlr] = tensorDenoiseSVD(DataT.Ysm, 0.98);
mlr

%% cut trials
rr = 10;
minTrial = min(Data.trialCount,[],2);
n_inds = minTrial >= rr;
n_inds = n_inds(1:size(Data.Y,1));

Ys = Data.Ys(n_inds,:,:,1:rr);

Yst = Ys(:,:,:);

%% denoise
% good neurons for DataLaraB: 4, 10, 27
n_inds = [4, 10, 27];
cplot = randperm(c, 2);
rplot = randperm(rr, rr);
Yst_ = tensorDenoiseSVD(Yst, mlr);
Yst_ = reshape(Yst_, size(Ys));
tensorSubplot(Yst_(n_inds,:,cplot,rplot),3,1)


