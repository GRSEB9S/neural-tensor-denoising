function tensorDenoiseDat2Plot(datpath)
% process dat files and make a plot
% make sure datpath includes '/' at the end

files = dir(datpath);
files = files(~[files.isdir]);
n = length(files);

load('mat-files/Data.mat');
load('mat-files/DataSim.mat');
load('mat-files/DataMatched.mat');

% initialize errPlot with nans/ []'s

for nn = 1:n
  load([datpath files(nn).name]);
  
  Yest = table_out.Properties.UserData.Yest;
  dataset = table_out.Properties.UserData.dataset;
  
  rng_iter = str2double(files(nn).name(5:7));
  trial = str2double(files(nn).name(9:11));
  method = str2double(files(nn).name(13));
  dataset = str2double(files(nn).name(15));

  switch dataset
    case 1
      Ygt = Data.Ysm_;
    case 2
      Ygt = DataSim.Ygt;
    case 3
      Ygt = DataMatched.Ysm_;
  end
  
  err = norm(Ygt(:) - Yest(:)).^2/Ygt(:).^2;
  
  errPlot{dataset}{rng_iter}(method+1,trial-1) = err;
    
end
end
