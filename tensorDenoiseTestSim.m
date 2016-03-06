%% tensor denoising test script using simulated data

% data types:
% simultaneous{c}(n,t,r) -- all r roughly equal and large
% simultaneous{c}(n,t,r) -- some conditions w/ low r, some with high r

% sequential{n}{c}(t,r)

% simultaneous denoising scenarios:
% [1] Treat data as (n,t,c). Avg across trials. Denoise filtered data with HOSVD.
% [2] Treat data as (n,t,R). Learn condition space.
% [2]b Treat data as (n,t,c). Learn core tensor. Project trials onto core tensor. 
% [3] c x c  structuured data fusion




%%
clear; clc; close all;

%% simulated data
[ Yout, Xsim, Bsim, s, f] = simulateTensorSpikes;
%% pick data type
Y = Yout.Y;
Yfilt = filterGauss(Y,10,2);
%% visualize
tensorSubplot(mean(Y,4,'omitnan'))
tensorSubplot(mean(Yfilt,4,'omitnan'))
tensorSubplot(Xsim)
















%% [1] - hosvd
Y1 = mean(Y,4,'omitnan');

options = [];
options.XMultiplier = 10;
options.Range = -1;
[sizeCore,L_ub] = mlrankest(Y1,options);

options = [];
options.Display = true;
options.AlgorithmOptions.TolFun = 1e-4;
options.AlgorithmOptions.TolX = 1e-4;
[Uhat, Shat, output] = lmlra(Y1, sizeCore, options);

[Uhat2,Shat2] = mlsvd(Y1,sizeCore);

Y1hat = lmlragen(Uhat, Shat);
Y1hat2= lmlragen(Uhat2, Shat2);

%% [1] - admm
[X, B] = tensorDenoiseADMM(Y, 'maxIter', 10, 'bf', 10, 'lf', 10, 'ef', 1, 'trueX', Xsim, 'trueB',Bsim);

%% [2] - n x t x R
Y2 = Yfilt(:,:,:);

%% set trank manually
trank = [10 6 6];

%% mlrankest (expensive)
trank = mlrankest(Y2);

%% 
[U,S,output] = lmlra(Y2,trank);
Y2est1 = reshape(tmprod(S,U,1:3), size(Yfilt));

[U,S,sv] = mlsvd(Y2,trank);
Y2est2 = reshape(tmprod(S,U,1:3), size(Yfilt));

%% plot a few trials
tensorSubplot(squeeze(Y2est1(:,:,1:3,:)))

%% [3] - structured data fusion. c x c across-data similarity

%% model struct
trank1 = [10 6 6];
trank2 = [];
model.variables.u = [];







































