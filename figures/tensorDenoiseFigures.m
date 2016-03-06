%% tensor denoising figures

clear; clc; close all;
%% load data from tensorDenoiseMakeData.m
load('TensorDenoising/2015/mat files/Data.mat');
load('TensorDenoising/2015/mat files/DataFull.mat');
load('TensorDenoising/2015/mat files/DataMatched.mat');
load('TensorDenoising/2015/mat files/DataSim.mat');

%% ------------------------------------------------------------------------
%% load previous sims
real = load('/Users/Jeff/Documents/MATLAB/TensorDenoising/2015/mat files/02-10-2016-real.mat');
simulated = load('/Users/Jeff/Documents/MATLAB/TensorDenoising/2015/mat files/02-09-2016-simulated.mat');
fakereal = load('/Users/Jeff/Documents/MATLAB/TensorDenoising/2015/mat files/02-10-2016-fakereal.mat');

%% figure 1 - example neurons
% good neurons. 97, 23, 37 (87, 10)...101?
figap = {...
         'dots',51,...
         'ps',[100 60],...
         'ybar',20,...
         'xbar',20,...
         'sort',30,...
        };
nn = 99;
[h1, iis] = figA(tensorDenoiseSVD(1000*Data.Ysm, 0.975), nn, figap{:});figPapersize
h2 = figA(1000*Ysm, nn, figap{:},'iis',iis);figPapersize

yl = [min(h1.Children.YLim(1), h2.Children.YLim(1)), max(h1.Children.YLim(2), h2.Children.YLim(2))];
h1.Children.YLim = yl;
h2.Children.YLim = yl;

%% figure 1 -- alternative. Show 2 example conditions / make sure summary is loaded
% neurons (after sortType = trialSort): 14 (note, 14 & 15 are identical?)
% conditions: 7 and 17.
figap = {...
         'dots',51,...
         'ps',[100 60],...
         'ybar',20,...
         'xbar',20,...
         'sort',30,...
        };
nn = 14;
cc = [7 17];
r1 = 3;
r2 = 20;
YDN = tensorDenoiseSVD(mean(Ys_(:,:,:,1:r1),4,'omitnan'), minrank(r1,:));
[h1, iis] = figA(1000*mean(Ys_(:,:,cc,1:r1),4,'omitnan'), nn, figap{:}); figPapersize
h2 = figA(1000*mean(Ys_(:,:,cc,1:r2),4,'omitnan'), nn, figap{:}); figPapersize
h3 = figA(1000*YDN(:,:,cc), nn, figap{:}); figPapersize

%% figure 2 - err vs trials / grid search / loocv /
% implement: bootstrap
options = [];
options.gridStep = 3;
options.minRank = [20 6 6];
options.maxRank = [80 20 20];
options.maxr = 10;
options.resample = 1;

Ys_ = tensorDenoiseStandardize(Data.Ys_);
Ytrue = mean(Ys_,4,'omitnan');

real(1).summary = tensorDenoiseGridSearchCV(Ys_, options );
for mm = 1:3
  options.mode = mm; real(mm+1).summary = tensorDenoiseUnfoldingCV(Ys_, options );
end
real(5).summary = tensorDenoiseFilterCV(Ys_, options );

%% plot results
tensorDenoiseSummaryPlot(real, Ytrue, Ys_);

%% figure 3 - err vs trials / grid search / loocv / simulated data
options = [];
options.gridStep = 3;
options.minRank = [20 6 6];
options.maxRank = [80 20 20];
options.maxr = 10;
options.resample = 1;

[Ytrue, mu, sig] = tensorDenoiseStandardize(DataSim.Ygt);
Ys_ = bsxfun(@plus, DataSim.Ys_, -mu);
Ys_ = bsxfun(@rdivide, Ys_, sig);

sim(1).summary = tensorDenoiseGridSearchCV( Ys_, options );
for mm = 1:3
  options.mode = mm; sim(mm+1).summary = tensorDenoiseUnfoldingCV( Ys_, options );
end
sim(5).summary = tensorDenoiseFilterCV( Ys_, options );


%% plot results
tensorDenoiseSummaryPlot(sim, Ytrue, Ys_, DataSim.size_core);

%% figure 3 - err vs trials / grid search / loocv / 'fake-real' data
options = [];
options.gridStep = 3;
options.minRank = [20 6 6];
options.maxRank = [80 20 20];
options.maxr = 10;
options.resample = 1;

Ys_ = tensorDenoiseStandardize(DataMatched.Ys_);
Ytrue = mean(Ys_,4,'omitnan');

matched(1).summary = tensorDenoiseGridSearchCV( Ys_, options );
for mm = 1:3
  options.mode = mm; matched(mm+1).summary = tensorDenoiseUnfoldingCV( Ys_, options );
end
matched(5).summary = tensorDenoiseFilterCV( Ys_, options );


%% plot results
tensorDenoiseSummaryPlot(matched, Ytrue, Ys_);




%% figure 4 -- missing data...
%% missing data
Ymiss = Ysm_;
m_inds = randi(maxn*c, 0.5*maxn*c, 1);
Ymiss = permute(Ymiss, [2 3 1]);
Ymiss(:,m_inds) = nan;
Ymiss = permute(Ymiss, [3 1 2]);

%% 

%%
options = [];
options.Display = true;
options.Algorithm = @lmlra_nls;
options.AlgorithmOptions.TolFun = 1e-8;
options.AlgorithmOptions.TolX = 1e-2;
options.AlgorithmOptions.Display = 1;
options.AlgorithmOptions.MaxIter = 50;

[U,S] = lmlra(Ymiss, [30 10 10], options);

Yrecon = lmlragen(U,S);


%% leave neuron out generalization
A1 = tensorDenoiseSVD(randn(10,10,10), [5 5 5]);
A2 = A1;
A2(randperm(numel(A2), round(0.7*numel(A2)))) = nan;
A2 = fmt(A2);

options = [];
options.Display = true;

[U,S] = lmlra(A2, [5 5 5], options);
A3 = lmlragen(U,S);

tensorSubplot(A1);
%tensorSubplot(A2);
tensorSubplot(A3);

options = [];
options.threshold = 0;
tensorDenoiseERR(A1,A3,options)

%% figure 2c - bonus analysis
Y_ = Ys(:,:,:,1:5);

Ysvd = tensorDenoiseSVD(mean(Y_,4,'omitnan'), 0.99);
% now add neuron and condition
maxr = 50;
maxn = 10;

Y__ = zeros(maxn,t,c,maxr);
newTrials = [];
for nn = 1:maxn
  for cc = 1:c
    newTrials{nn,cc} = randi([6 trialCount(nn,cc)],maxr,1);
    Y__(nn,:,cc,:) = Ys(nn,:,cc,newTrials{nn,cc});
  end
end
Ypsth = cat(4, Y_(1:10,:,:,:), Y__);

%% figure 2c - calculate error and plot
err = [];
for rr = 5:size(Y__,4)
  err{rr} = tensorDenoiseERR(Ysm(1:10,:,:,:), mean(Ypsth(:,:,:,1:rr), 4, 'omitnan'), 1, 2);
end
err = [err{:}];
errRef = tensorDenoiseERR(Ysm(1:10,:,:,:), Ysvd(1:10,:,:,:), 1, 2);

figure; hold all;
nn = 6;

plot(5:maxr, err(nn,:),'b')
plot([5 maxr], [errRef(nn) errRef(nn)]);

%%






