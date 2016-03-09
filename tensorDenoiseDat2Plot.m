function tensorDenoiseDat2Plot(datpath)
% process dat files and make a plot
% make sure datpath includes '/' at the end

%%%% REDO ...

files = dir(datpath);
files = files(~[files.isdir]);
n = length(files);

load('mat-files/Data.mat');
load('mat-files/DataSim.mat');
load('mat-files/DataMatched.mat');

% initialize errPlot with nans/ []'s

for nn = 1:n
  load([datpath files(nn).name]);
  
  Yest = table_out.Yest;
  %dataset = table_out.dataset;
  
  rng_iter = str2double(files(nn).name(5:7));
  trial = str2double(files(nn).name(9:11));
  method = str2double(files(nn).name(13));
  dataset = str2double(files(nn).name(15));

  switch dataset
    case 1
      Ygt = Data.Ysm_;
      [~,mu,sig] = tensorDenoiseStandardize(Data.Ys_);
      Ygt = bsxfun(@plus, Ygt, -mu);
      Ygt = bsxfun(@rdivide, Ygt, sig);
    case 2
      Ygt = DataSim.Ygt;
      [~,mu,sig] = tensorDenoiseStandardize(DataSim.Ygt);
      Ygt = bsxfun(@plus, Ygt, -mu);
      Ygt = bsxfun(@rdivide, Ygt, sig);
    case 3
      Ygt = DataMatched.Ysm_;
      [~,mu,sig] = tensorDenoiseStandardize(DataMatched.Ys_);
      Ygt = bsxfun(@plus, Ygt, -mu);
      Ygt = bsxfun(@rdivide, Ygt, sig);
  end
  
  err = norm(Ygt(:) - Yest(:)).^2./norm(Ygt(:)).^2;
  
  errPlot{dataset}(trial-1,method+1,rng_iter+1) = err;
  minRank{dataset}{method+1}(trial-1,:,rng_iter+1) = table_out.minrank;
    
end

figure; hold all
xaxis = 2:size(errPlot{1},1)+1;
plot(xaxis,errPlot{1}(:,:))
legend('tensor','neuron','time','condition');
xlabel('trials (per neuron per condition)');
ylabel('relative error');
title('main analysis | relative error plot')

figure; hold all
plot(xaxis,minRank{1}{1}(:,:))


figure; hold all
xaxis = 2:size(errPlot{1}{1},2)+1;
plot(xaxis,errPlot{1}{1}')
legend('tensor','neuron','time','condition');
xlabel('trials (per neuron per condition)');
ylabel('relative error');
title('main analysis | relative error plot')

end
