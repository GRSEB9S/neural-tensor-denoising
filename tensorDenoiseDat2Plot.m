function [errPlot, minRankPlot] = tensorDenoiseDat2Plot(datindex)
% process dat files and make a plot
% make sure datpath includes '/' at the end

%%
if nargin == 0 
  datindex = 1;
end

if isdir('/ifs/')
  cd /ifs/scratch/zmbbi/la_lab/jss2219/
  addpath(genpath('/ifs/scratch/zmbbi/la_lab/jss2219/'))
elseif isdir('/vega/')
  cd /vega/zmbbi/users/jss2219/TensorDenoising/
  addpath(genpath('/vega/zmbbi/users/jss2219/TensorDenoising/'))
end

files_r = dir('dat-file*');
datpath = [files_r(datindex).name '/'];

files = dir(datpath);
files = files(~[files.isdir]);
n = length(files);
load([datpath files(1).name]);
load(summary.options.refData);

if summary.options.simrun
  for dd = 1:summary.options.n_data
    D(dd).Data = tensorDenoiseMakeDataSim(Data, summary.options.rnks(dd,:));
  end
end

%% get errPlot and minRankPlot 
% initialize errPlot with nans/ []'s
errPlot = cell(summary.options.n_data, 1);
errPlot(:) = {nan(summary.options.n_trialcount, summary.options.n_method, summary.options.n_iter)};

minRankPlot = cell(summary.options.n_data, summary.options.n_method);
minRankPlot(:) = {nan(summary.options.n_trialcount, size(summary.minrank,2), summary.options.n_iter)};

for nn = 1:n
  load([datpath files(nn).name]);
  
  r1 = summary.options.r1;
  rng_iter = summary.options.rng_iter;
  method = summary.options.method;
  minrank = summary.minrank;
  dataset = summary.options.dataset;
  %minrank = summary.minrank_p(20,:);

  if summary.options.simrun
    Ygt = D(dataset).Data.Ygt;
    Y = resampleTrials(D(dataset).Data.Ys, 1, r1+1000*rng_iter);    
  else
    Ygt = Data.Ysm;
    Y = resampleTrials(Data.Ys, 1, r1+1000*rng_iter);
  end
  Yest = tensorDenoiseSVD(mean(Y(:,:,:,1:r1),4,'omitnan'), minrank);
  Yest = Yest.*(Yest>0);
  
  err = norm(Ygt(:) - Yest(:)).^2./norm(Ygt(:)).^2;
  
  errPlot{dataset}(r1-1,method+1,rng_iter+1) = err;
  minRankPlot{dataset, method+1}(r1-1,:,rng_iter+1) = summary.minrank;
    
end

%% average method
maxrng_iter = summary.options.n_iter;
maxtrial = summary.options.n_trialcount;
maxmethod = summary.options.n_method;
maxdataset = summary.options.n_data;

for dataset = 1:maxdataset
  for r1 = 2:maxtrial+1 % 
    for rng_iter = 0:maxrng_iter
        Ygt = Data.Ysm;
        Y = resampleTrials(Data.Ys, 1, r1+1000*rng_iter); % probably can put this outside of the loop
        Yest = mean(Y(:,:,:,1:r1),4,'omitnan');
        err = norm(Ygt(:)-Yest(:)).^2./norm(Ygt(:)).^2;
        errPlot{dataset}(r1-1,maxmethod+1,rng_iter+1) = err;
    end
  end
end

%% save
save([datpath 'plotfiles'], 'errPlot','minRankPlot');

end
