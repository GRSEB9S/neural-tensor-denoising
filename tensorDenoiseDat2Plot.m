function tensorDenoiseDat2Plot(datpath)
% process dat files and make a plot
% make sure datpath includes '/' at the end

files = dir(datpath);
files = files(~[files.isdir]);
n = length(files);

for nn = 1:n
  rng_iter = str2num(files(nn).name(?));
  trial = str2num(files(nn).name(?));
  method = str2num(files(nn).name(?));
  dataset = str2num(files(nn).name(?));
  
  dat = readtable([datpath files(nn).name]);
  
  minrank = dat{1,{'minrank_1', 'minrank_2', 'minrank_3'}};
  err = dat{:,'err'};
  ranks = dat{:,{'ranks_1','ranks_2','ranks_3'}};
  minidx = find(ismember(ranks,minrank,'rows')==1);
  
  errPlot(dataset).method(method+1).err(trial-1,1) = err(minidx); % use cell instead? pad/initialize with nans in case rows/cols don't match up
  
end

end
