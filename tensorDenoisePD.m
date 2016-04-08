%% compute preferred direction
% take in (n,t,c) tensor.
% for each (n,t), compute the PD.
% possibly consider longer time windows, e.g. sliding time bin averaging.

%% init
clear;clc;close all;

%%
load('/Users/Jeff/Documents/MATLAB/TensorDenoising/mat-files/DataLaraACat')
% for DataLaraCat:
% target on: 400ms (40)
% end targlock: 900ms
% begin movelock: 901ms
% move: 1400ms (140)
%

r1a = 3;
r1b = 20;

[n,t,c,r] = size(Data.Y);
[~,varSort] = sort(var(Data.Ysm(:,:),[],2),'descend');

%% early trials vs late trials
Data2 = Data;
for nn = 1:n
  for cc = 1:c
    Data2.Y(nn,:,cc,1:Data.trialCount(nn,cc)) = Data.Y(nn,:,cc,Data.trialCount(nn,cc):-1:1);
    Data2.Ys(nn,:,cc,1:Data.trialCount(nn,cc)) = Data.Ys(nn,:,cc,Data.trialCount(nn,cc):-1:1);
    Data2.Ysm(nn,:,cc,1:Data.trialCount(nn,cc)) = Data.Ys(nn,:,cc,Data.trialCount(nn,cc):-1:1);
  end
end

%% Load target data
% NA or NB
load('/Users/Jeff/Documents/MATLAB/Datasets/Motor/Lara/NA.mat');

%% get targ values;
% note: each N(1), N(2), etc, has a different x,y, coordinate for the
% targets, but seem to be essentially at the same angle. so as long as we
% do arctan we should be fine(?)
targ = zeros(2,c);
for cc = 1:c
  targ(1,cc) = N(1).condTable(cc,5).TARGINFO.x/150; % put in sensible range /150.
  targ(2,cc) = N(1).condTable(cc,5).TARGINFO.y/150;
end
targ = cat(1,targ, ones(1,c));

%% get denoise summary
options = struct;
options.gridStep = 4;
options.minRank = [40 6 6];
options.maxRank = [80 30 22];
options.resample = 0;
options.verbose = 1;

options.rng_iter = 0;
options.r1 = r1a;
options.method = 0;

% summarya = tensorDenoiseGridSearchCV(Data.Ys, options);

%% r1b trials
% options.r1 = r1b;
% summaryb = tensorDenoiseGridSearchCV(Data.Ys, options);

%% r1a trials / late
% options.r1 = r1a;
% summaryc = tensorDenoiseGridSearchCV(Data2.Ys, options);

%% denoise
% Yhata = tensorDenoiseSVD(Data.Ysm(:,:,:,1:r1a), summarya.minrank);
% Yhatb = tensorDenoiseSVD(Data.Ysm(:,:,:,1:r1b), summaryb.minrank);
% Yhatc = tensorDenoiseSVD(Data2.Ysm(:,:,:,1:r1a), summaryc.minrank);

Yhat{1} = filterGauss(mean(Data.Ys(:,:,:,1:r1a),4,'omitnan'),3,2);
Yhat{2} = filterGauss(mean(Data.Ys(:,:,:,1:r1b),4,'omitnan'),3,2);
Yhat{3} = filterGauss(mean(Data2.Ys(:,:,:,1:r1a),4,'omitnan'),3,2);

Yhat{4} = tensorDenoiseSVD(mean(Data.Ys(:,:,:,1:r1a),4,'omitnan'), [20 10 5]);
Yhat{5} = tensorDenoiseSVD(mean(Data.Ys(:,:,:,1:r1b),4,'omitnan'), [20 10 5]);
Yhat{6} = tensorDenoiseSVD(mean(Data2.Ys(:,:,:,1:r1a),4,'omitnan'), [20 10 5]);

Yhat{7} = filterGauss(Data.Ysm, 3, 2);

%% get PD
for yy = 1:7
  pd(:,:,yy) = tensorDenoiseGetPD(Yhat{yy}, targ);
end

%% plot
% good inds for DataLaraA:
n_inds = [2 5 6 7 8 9 11 13 16 17 18 20 21 23];
%n_inds = 1:20;

cplot = [1 4 7];
h = tensorSubplot(pd(n_inds,:,cplot),4,2,'csort',0);
xlabel('time');
ylabel('PD (degrees)');
for aa = 1:length(h.Children)
  h.Children(aa).XTick = [40 140];
  h.Children(aa).XTickLabel = {'on','move'};
end

% target on: 400ms (40)
% end targlock: 900ms
% begin movelock: 901ms
% move: 1400ms (140)

legend({'avg+filter, first 3 trials','tensor denoise first 3 trials','avg+filter, all trials (~ ground truth)'});

%legend({'average (10ms filter)', 'average (+50ms filter)', 'tensor denoised (CV rank)', 'tensor denoised (oversmoothed)'});
h = tensorSubplot(Yhat{7}(n_inds,:,:),4,2);
for aa = 1:length(h.Children)
  h.Children(aa).XTick = [40 140];
  h.Children(aa).XTickLabel = {'on','move'};
end

%%

