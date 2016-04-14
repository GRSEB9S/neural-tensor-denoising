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

r1a = 5;
r1b = 20;

%% cut neurons to high-r
[n,t,c,r] = size(Data.Y);
[~,varSort] = sort(var(Data.Ysm(:,:),[],2),'descend');

minTrial = min(Data.trialCount,[],2);
minTrial = minTrial(1:size(Data.Ys,1));
n_inds = find(minTrial>=10);

Data.trialCount = Data.trialCount(n_inds,:);
r = max(Data.trialCount(:));
n = length(n_inds);

Data.Y = Data.Y(n_inds,:,:,1:r);
Data.Ys = Data.Ys(n_inds,:,:,1:r);
Data.Ysm = Data.Ysm(n_inds,:,:);

%% resampled trials
Datax = Data;
Datax.Y = resampleTrials(Data.Y, 0, 1);
Datax.Ys = resampleTrials(Data.Ys, 0, 1);
Datax.Ysm = resampleTrials(Data.Ysm, 0, 1);

%% early trials vs late trials
Data2 = Data;
Datax2 = Datax;
for nn = 1:n
  for cc = 1:c
    Data2.Y(nn,:,cc,1:Data.trialCount(nn,cc)) = Data.Y(nn,:,cc,Data.trialCount(nn,cc):-1:1);
    Data2.Ys(nn,:,cc,1:Data.trialCount(nn,cc)) = Data.Ys(nn,:,cc,Data.trialCount(nn,cc):-1:1);
    Data2.Ysm(nn,:,cc,1:Data.trialCount(nn,cc)) = Data.Ys(nn,:,cc,Data.trialCount(nn,cc):-1:1);
    
    Datax2.Y(nn,:,cc,1:Data.trialCount(nn,cc)) = Datax.Y(nn,:,cc,Data.trialCount(nn,cc):-1:1);
    Datax2.Ys(nn,:,cc,1:Data.trialCount(nn,cc)) = Datax.Ys(nn,:,cc,Data.trialCount(nn,cc):-1:1);
    Datax2.Ysm(nn,:,cc,1:Data.trialCount(nn,cc)) = Datax.Ys(nn,:,cc,Data.trialCount(nn,cc):-1:1);
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
mlr = [20 10 10];
filt = 3;

% Yhata = tensorDenoiseSVD(Data.Ysm(:,:,:,1:r1a), summarya.minrank);
% Yhatb = tensorDenoiseSVD(Data.Ysm(:,:,:,1:r1b), summaryb.minrank);
% Yhatc = tensorDenoiseSVD(Data2.Ysm(:,:,:,1:r1a), summaryc.minrank);

Yhat{1} = filterGauss(mean(Data.Ys,4,'omitnan'),filt,2);

Yhat{2} = filterGauss(mean(Data.Ys(:,:,:,1:r1a),4,'omitnan'),filt,2);
Yhat{3} = filterGauss(mean(Data2.Ys(:,:,:,1:r1a),4,'omitnan'),filt,2);

Yhat{4} = filterGauss(mean(Datax.Ys(:,:,:,1:r1a),4,'omitnan'),filt,2);
Yhat{5} = filterGauss(mean(Datax2.Ys(:,:,:,1:r1a),4,'omitnan'),filt,2);

Yhat{6} = tensorDenoiseSVD(mean(Data.Ys(:,:,:,1:r1a),4,'omitnan'), mlr);
Yhat{7} = tensorDenoiseSVD(mean(Data2.Ys(:,:,:,1:r1a),4,'omitnan'), mlr);

Yhat{8} = tensorDenoiseSVD(mean(Datax.Ys(:,:,:,1:r1a),4,'omitnan'), mlr);
Yhat{9} = tensorDenoiseSVD(mean(Datax2.Ys(:,:,:,1:r1a),4,'omitnan'), mlr);

datasets = {...
            'all trials',...
            'avg, first block',...
            'avg, last block',...
            'avg, rand first block',...
            'avg, rand last block',...
            'tensor, first block',...
            'tensor, last block',...
            'tensor, rand first block',...
            'tensor, rand last block',...
            };

%% get tuning strength
ts = [];
for yy = 1:length(Yhat)
  ts(:,:,yy) = var(Yhat{yy}, 0, 3);
end
%% plot tuning
dplot = [1 6 7];
n_inds = [2 5 6 7 8 9 11 13 16 17 18 20 21 23];
tensorSubplot(ts(n_inds, :, dplot), 4,2,'csort',0);

%% get PD
pd = [];
for yy = 1:length(Yhat)
  pd(:,:,yy) = tensorDenoiseGetPD(Yhat{yy}, targ);
end

%% plot
% good inds for DataLaraA:
%n_inds = [2 5 6 7 8 9 11 13 16 17 18 20 21 23];
%n_inds = 1:20;

dplot = [1 4 7];
tensorSubplot(ts(n_inds, :, dplot), 4,2,'csort',0);
h = tensorSubplot(pd(n_inds,:,dplot),4,2,'csort',0);
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

legend(datasets{dplot});

h = tensorSubplot(Yhat{1}(n_inds,:,:),4,2);
for aa = 1:length(h.Children)
  h.Children(aa).XTick = [40 140];
  h.Children(aa).XTickLabel = {'on','move'};
end

%% some analysis... find strongest tuned part of psth
tuneTime = [];
timeRange = 120:180;
yref = 1;
for nn = 1:n
  for yy = 1:length(Yhat)
    [~,tuneTime(nn,yy)] = max(ts(nn,timeRange,yy));
  end
end
tuneTime = tuneTime + timeRange(1) - 1;

%% calculate pd diff
yref = ones(length(Yhat),1);
yref = [1 3 3 5 5 7 7 9 9];

for nn = 1:n
  for yy = 1:length(Yhat)
    pddiff_(nn,yy) = pd(nn,tuneTime(nn,yref(yy)),yy) - pd(nn,tuneTime(nn,yref(yy)),yref(yy));
  end
end
pddiff = min(abs(pddiff_), 360 - abs(pddiff_));

%% get varsort
for nn = 1:n
  tssort(nn) = ts(nn,tuneTime(nn,yref(1)),yref(1));
end
[~,varsort] = sort(tssort,'descend');
varsort = varsort(1:50);
meanpddiff = mean(pddiff);

%% plot
dplot = [2 4 6 8];
h=figure;
plot(pddiff(varsort,dplot));
hold all
h.Children(1).ColorOrderIndex = 1;
for dd = 1:length(dplot)
  plot([1 50], [meanpddiff(dplot(dd)) meanpddiff(dplot(dd))])
end
legend(datasets{dplot});























