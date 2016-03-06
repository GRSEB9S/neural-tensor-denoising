%% tensor denoising tests on real data
% test tensor denoising on M1, single trial, and olfaction data

% step 1: basic hosvd on all data
% step 2: basic hosvd for single trials
% step 3: admm on all data
% step 4: admm on single trials
%
% sdf steps








%%
clear; clc; close all;

%% paths
addpath('Datasets/Axel (Olfaction)/')
addpath('#Packages/Seely - Tensor Denoising/');
addpath('TensorDenoising/2015/');

%% load data
Shenoy = preRC2('Datasets/Motor/Shenoy/RC,2009-09-18,1-2,good-ss', 10, 10);
Lara = preS2('Datasets/Motor/Lara-20141105/Sstructs/', 10, 10);
%Axel = axelFun('binWidth',100,'filterWidth',150);

%% select shenoy alignment
Shenoy2 = Shenoy.aligned(3);
Shenoy2.cond = Shenoy.cond;

%% format data
ShenoyFormatted = compressTrials2(Shenoy2);
LaraFormatted = compressTrials2(Lara);
%AxelFormatted = compressTrials2(Axel);

%% select data
Y = ShenoyFormatted.A;
Ys = ShenoyFormatted.As;
[~,sortn] = sort(sum(Y(:,:)','omitnan'),'descend');
Y = Y(sortn,:,:,:);
Ys = Ys(sortn,:,:,:);





%% [1] - ml rank - tensorlab
Y1 = mean(Ys,4,'omitnan');

options = [];
options.XMultiplier = 10;
options.Range = -1;
[sizeCore,L_ub] = mlrankest(Y1,options);

sizeCore = L_ub(650,1:3);

options = [];
options.Display = true;
options.AlgorithmOptions.TolFun = 1e-4;
options.AlgorithmOptions.TolX = 1e-4;
[Uhat, Shat, output] = lmlra(Y1, sizeCore, options);

[Uhat2,Shat2] = mlsvd(Y1,sizeCore);

Y1hat = lmlragen(Uhat, Shat);
Y1hat2= lmlragen(Uhat2, Shat2);

%% plot // frob(T - T_hat) 
  semilogy(0:output.Algorithm.iterations,sqrt(2*output.Algorithm.fval));
  xlabel('iteration');
  ylabel('frob(lmlrares(T,U))');
  grid on;

%% compute the residual // compute the error
  res = lmlrares(Y1, Uhat, Shat);
  relerr = frob(lmlrares(Y1, Uhat, Shat))/frob(Y1);
  

%% [2]
Y2 = Ys(:,:,:);

[U,S] = lmlra(Y2,sizeCore);
Y2est = reshape(tmprod(S,U,1:3), size(Y2));

















%% downselect neurons
%[~,sortn] = sort(sum(ShenoyFormatted.A(:,:)','omitnan'),'descend');
%ShenoyDownselected = ShenoyFormatted.A(sortn(1:50),:,:,:);

%% denoise data / formatted / hosvd
r = 40;
ShenoyFormattedHOSVD = tensorDenoise(nanmean(ShenoyFormatted.As,4), r);

%% denoise data / unformatted / hosvd
r = 80;
pre = 30;
post = 30;
for ii = 1:4
   timeZero = Shenoy.aligned(ii).timeZero;
   t1 = timeZero - pre;
   t2 = timeZero + post;
   [ShenoyHOSVD{ii}.As, ShenoyHOSVD{ii}.hosvd_dims] = tensorDenoise(Shenoy.aligned(ii).As(:,t1:t2,:), r);
   ShenoyHOSVD{ii}.A = Shenoy.aligned(ii).A(:,t1:t2,:);
   ShenoyHOSVD{ii}.cond = Shenoy.cond;
end

%% reformat ShenoyHOSVD to 4th order tensor
for ii = 1:4
   timeZero = Shenoy.aligned(ii).timeZero;
   t1 = timeZero - pre;
   t2 = timeZero + post;
   alignedStruct = Shenoy.aligned(ii);
   alignedStruct.As = ShenoyHOSVD{ii};
   alignedStruct.cond = Shenoy.cond;
   ShenoyFormattedHOSVD{ii} = compressTrials2(alignedStruct);
end
%% plot reformatted data
[~,sortn] = sort(sum(Shenoy.aligned(4).A(:,:)', 'omitnan'), 'descend');
tensorSubplot(mean(ShenoyFormattedHOSVD{4}.As)
for cc = [1 11 21 31 41]
   [~,sortn] = sort(sum(ShenoyFormattedHOSVD{4}.A(:,:)','omitnan'),'descend');
   tensorSubplot(squeeze(ShenoyFormattedHOSVD{4}.As(sortn,:,cc,:)))
end



%% ADMM denoise / simulated data
[X, B] = tensorDenoiseADMM(Ysim,'maxIter',10,'bf',3,'lf',1,'ef',1,'trueX',Xsim,'trueB',Bsim);
%% plot shit
data{1} = X;
data{2} = Xsim;
tensorComparePlot(data,10);

%%
tensorSubplot(X);
tensorSubplot(Xsim);
tensorSubplot(mean(Ysim,4,'omitnan'))







































