function tensorDenoiseCluster(index)
%% generate figure 2 data. for hpc cluster.

options = struct;

%% set options
%options.refData = './mat-files/DataShenoy.mat';
options.refData = './mat-files/DataLaraA.mat';
options.datpath = './dat-file-Lara-8-17-Asim/';
options.simrun = 1;

% basic gridsearchCV options
options.gridStep = 4;
options.minRank = [40 6 6];
options.maxRank = [80 30 22];
options.resample = 1;
options.verbose = 1;

%% parameters for iteration
if options.simrun
  % sim data ranks
  rnk_grid{1} = [30];
  rnk_grid{2} = [15];
  rnk_grid{3} = [10];
  rnks = cell(1,3);
  [rnks{:}] = ndgrid(rnk_grid{:});
  options.rnks = cell2mat(cellfun(@(i)i(:),rnks,'UniformOutput',false));
  options.n_data = size(options.rnks,1);
else
  options.n_data = 1;
end

options.n_method = 5;
options.n_trialcount = 20;
options.n_iter = 5;
n_vec = [options.n_iter, options.n_trialcount, options.n_method, options.n_data];
disp(prod(n_vec)) % display index range. 

%% paths 
if isdir('/ifs/')
  cd /ifs/scratch/zmbbi/la_lab/jss2219/
  mkdir(options.datpath);
  addpath(genpath('/ifs/scratch/zmbbi/la_lab/jss2219/'))
elseif isdir('/vega/')
  cd /vega/zmbbi/users/jss2219/TensorDenoising/
  mkdir(options.datpath);
  addpath(genpath('/vega/zmbbi/users/jss2219/TensorDenoising/'))
else
  cd /Users/jeff/Documents/MATLAB/TensorDenoising/
end

load(options.refData);

%% array of parameters
in = arrayfun(@(i) 1:n_vec(i), 1:length(n_vec),'UniformOutput',false);
out = cell(1,length(n_vec));
[out{:}] = ndgrid(in{:});
out = cell2mat(cellfun(@(i)i(:),out,'UniformOutput',false));

%% which one to perform
options.rng_iter = out(index,1)-1; % start with 0 
options.r1 = out(index,2)+1; % important: start with 2, not 1.
options.method = out(index,3)-1; % important: start with 0;
options.dataset = out(index,4);  

%% Run 
% make custom dataset if performing sim data
if options.simrun
  options.simrank = options.rnks(options.dataset, :);
  Data = tensorDenoiseMakeDataSim(Data, options.simrank);
end

% main grid search.
summary = tensorDenoiseGridSearchCV(Data.Ys, options);

% save 
filename = [options.datpath 'dat-' num2str(options.rng_iter, '%03d') '-' num2str(options.r1, '%03d') '-' num2str(options.method) '-' num2str(options.dataset, '%03d')];
save(filename, 'summary'); 

if options.verbose
  disp('Done')
end

quit force

end
