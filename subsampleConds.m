%% subsample conds

load('./mat-files/DataLaraA');

%%
% Data 1 - 4 trials / cond no T smooth
% Data 2 - 4 trials / cond 4 conds, 100 neurons, TS
% Data 3 - 2 trials / cond, 8 conds, 100 neurons, TS

%%
conds = randperm(Data.size(3), 8);
conds1 = conds(1:4);

X1 = Data.Ys(:,:,conds,1:4);
X2 = Data.Ys(:,:,conds1,1:4);
X3 = Data.Ys(:,:,conds,1:3);

Ygt = Data.Ysm(:,:,:);

%% tensor denoising
%% gs options
options = [];
options.gridStep = 0.01;
options.minRank = 0.7;
options.maxRank = 0.99;
options.resample = 1;
options.verbose = 1;

options.rng_iter = 0; % start with 0 
options.method = 4; % important: start with 0;
options.dataset = 0;  

%%
options.r1 = 4; % important: start with 2, not 1.
summary2 = tensorDenoiseGridSearchCV(X2, options);
options.r1 = 2; % important: start with 2, not 1.
summary3 = tensorDenoiseGridSearchCV(X3, options);


%%
th = 0.95;
X1_ = mean(X1,4);
X2_ = tensorDenoiseSVD(mean(X2,4), summary2.minrank);
X3_ = tensorDenoiseSVD(mean(X3,4), summary3.minrank);

%%
options = [];
options.perNeuron = 0;
options.errMeasure = 1;
options.threshold = 1;
options.relative = 1;

%%
nn = 1:80;
cc = 1:4;
err1 = tensorDenoiseERR(Ygt(nn,:,cc), X1_(nn,:,cc), options);
err2 = tensorDenoiseERR(Ygt(nn,:,cc), X2_(nn,:,cc), options);
err3 = tensorDenoiseERR(Ygt(nn,:,cc), X3_(nn,:,cc), options);

%%
bar([err1,err2,err3])
ylabel('rel err');


