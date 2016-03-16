function tensorDenoiseCluster(index)
% generate figure 2 data. for hpc cluster.

if isdir('/ifs/')
  cd /ifs/scratch/zmbbi/la_lab/jss2219/
else
  cd /Users/Jeff/Documents/MATLAB/TensorDenoising/
end
datpath = 'dat-file-3-16-16/';
mkdir(datpath);
addpath(genpath('/ifs/scratch/zmbbi/la_lab/jss2219/'))

%%
n_data = 1;
n_methods = 4;
n_trialcount = 10;
n_iter = 5;

n_vec = [n_iter, n_trialcount, n_methods, n_data];

disp(prod(n_vec))

%%
in = arrayfun(@(i) 1:n_vec(i), 1:length(n_vec),'UniformOutput',false);
out = cell(1,length(n_vec));
[out{:}] = ndgrid(in{:});
out = cell2mat(cellfun(@(i)i(:),out,'UniformOutput',false));

%% which one to perform
options = [];

% iteration options
options.rng_iter = out(index,1)-1; % start with 0 
options.r1 = out(index,2)+1; % important: start with 2, not 1.
options.method = out(index,3)-1; % important: start with 0;
options.dataset = out(index,4);

% basic options
options.gridStep = 4;
options.minRank = [20 4 4];
options.maxRank = [100 50 50];
options.resample = 1;
options.verbose = 1;


%% load data
load ('mat-files/DataShenoy.mat');

  %% comment out - for simulated data of different ranks
%   Data = tensorDenoiseMakeDataSim(Data, size_core);
%   rnk_n = [10 60];
%   rnk_t = [5 15];
%   rnk_c = [10 30];
  %% // comment out

summary = tensorDenoiseGridSearchCV(Data.Ys, options);
filename = [datpath 'dat-' num2str(options.rng_iter, '%03d') '-' num2str(options.r1, '%03d') '-' num2str(options.method) '-' num2str(options.dataset)];
save(filename, 'summary'); 

%%

quit force

end
