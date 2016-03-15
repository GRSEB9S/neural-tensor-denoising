function tensorDenoiseDat2Plot(datpath)
% process dat files and make a plot
% make sure datpath includes '/' at the end


files = dir(datpath);
files = files(~[files.isdir]);
n = length(files);

load('mat-files/Data_maze.mat');
%load('mat-files/DataSim.mat');
%load('mat-files/DataMatched.mat');

% initialize errPlot with nans/ []'s

% replace 'trial' with 'r1' below.
for nn = 1:n
  load([datpath files(nn).name]);
  
  r1 = table_out.options.r1;
  rng_iter = table_out.options.rng_iter;
  minrank = table_out.minrank;
  minrank = table_out.minrank_p(20,:);
  
  %dataset = table_out.dataset;
  
  rng_iter = str2double(files(nn).name(5:7));
  trial = str2double(files(nn).name(9:11));
  method = str2double(files(nn).name(13));
  dataset = str2double(files(nn).name(15));

  switch dataset
    case 1
      Ygt = Data.Ysm_;
      [Y,mu,sig] = tensorDenoiseStandardize(Data.Ys_);
      Ygt = bsxfun(@plus, Ygt, -mu);
      Ygt = bsxfun(@rdivide, Ygt, sig);
      Y = resampleTrials(Y, 1, rng_iter);
      Yest = tensorDenoiseSVD(mean(Y(:,:,:,1:r1),4,'omitnan'), minrank);
    case 2
      [Ygt,mu,sig] = tensorDenoiseStandardize(DataSim.Ygt);
      Y = bsxfun(@plus, DataSim.Ys_, -mu);
      Y = bsxfun(@rdivide, Y, sig);
      Y = resampleTrials(Y, 1, rng_iter);
      Yest = tensorDenoiseSVD(mean(Y(:,:,:,1:r1),4,'omitnan'), minrank);
    case 3
      Ygt = DataMatched.Ysm_;
      [Y,mu,sig] = tensorDenoiseStandardize(DataMatched.Ys_);
      Ygt = bsxfun(@plus, Ygt, -mu);
      Ygt = bsxfun(@rdivide, Ygt, sig);
      Y = resampleTrials(Y, 1, rng_iter);
      Yest = tensorDenoiseSVD(mean(Y(:,:,:,1:r1),4,'omitnan'), minrank);
  end
  
  err = norm(Ygt(:) - Yest(:)).^2./norm(Ygt(:)).^2;
  
  errPlot{dataset}(trial-1,method+1,rng_iter+1) = err;
  minRankPlot{dataset}{method+1}(trial-1,:,rng_iter+1) = table_out.minrank;
    
end

%% get average method
allfiles = vertcat(files(:).name);
maxrng_iter = max(str2num(allfiles(:,5:7)));
maxtrial = max(str2num(allfiles(:,9:11)));
maxmethod = max(str2num(allfiles(:,13)));
maxdataset = max(str2num(allfiles(:,15)));

for r1 = 2:maxtrial
  for rng_iter = 0:maxrng_iter
      Ygt = Data.Ysm_;
      [Y,mu,sig] = tensorDenoiseStandardize(Data.Ys_);
      Ygt = bsxfun(@plus, Ygt, -mu);
      Ygt = bsxfun(@rdivide, Ygt, sig);
      Y = resampleTrials(Y, 1, rng_iter);
      Yest = mean(Y(:,:,:,1:r1),4,'omitnan');
      err = norm(Ygt(:)-Yest(:)).^2./norm(Ygt(:)).^2;
      errPlot{1}(r1-1,5,rng_iter+1) = err;
      
%       [Ygt,mu,sig] = tensorDenoiseStandardize(DataSim.Ygt);
%       Y = bsxfun(@plus, DataSim.Ys_, -mu);
%       Y = bsxfun(@rdivide, Y, sig);
%       Y = resampleTrials(Y, 1, rng_iter);
%       Yest = mean(Y(:,:,:,1:r1),4,'omitnan');
%       err = norm(Ygt(:)-Yest(:)).^2./norm(Ygt(:)).^2;
%       errPlot{2}(r1-1,5,rng_iter+1) = err;
      
%       Ygt = DataMatched.Ysm_;
%       [Y,mu,sig] = tensorDenoiseStandardize(DataMatched.Ys_);
%       Ygt = bsxfun(@plus, Ygt, -mu);
%       Ygt = bsxfun(@rdivide, Ygt, sig);
%       Y = resampleTrials(Y, r1, rng_iter);
%       Yest = mean(Y(:,:,:,1:r1),4,'omitnan');
%       err = norm(Ygt(:)-Yest(:)).^2./norm(Ygt(:)).^2;
%       errPlot{3}(trial-1,5,rng_iter+1) = err;
  end
end

%% plot
figure; hold all
xaxis = 2:size(errPlot{1},1)+1;
plot(xaxis,mean(errPlot{1},3))
legend('tensor','neuron','time','condition','average');
xlabel('trials (per neuron per condition)');
ylabel('relative error');
title('main analysis | relative error vs trial')

figure; hold all
plot(xaxis,mean(minRankPlot{1}{1},3))
legend('neuron','time','condition');
xlabel('trials (per neuron per condition)');
ylabel('rank');
title('main analysis | optimal rank vs trial')

end
