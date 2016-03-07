function tensorDenoiseCluster(index)
% generate figure 2 data. for hpc cluster.

n_data = 3;
n_methods = 4;
n_trialcount = 20;
n_vec = [n_trialcount, n_methods, n_data];

in = arrayfun(@(i) 1:n_vec(i), 1:length(n_vec),'UniformOutput',false);
out = cell(1,length(n_vec));
[out{:}] = ndgrid(in{:});
out = cell2mat(cellfun(@(i)i(:),out,'UniformOutput',false));

%% which one to perform
r1 = out(index,1)+1; % important: start with 2, not 1.
method = out(index,2)-1; % important: start with 0;
dataset = out(index,3);

%% load data from tensorDenoiseMakeData.m
switch dataset
  case 1
    load('/mat-files/Data.mat');
    Ys_ = tensorDenoiseStandardize(Data.Ys_);
    tensorDenoiseWrapper(Ys_, r1, method, dataset);
  case 2
    load('/mat-files/DataSim.mat');
    [~, mu, sig] = tensorDenoiseStandardize(DataSim.Ygt);
    Ys_ = bsxfun(@plus, DataSim.Ys_, -mu);
    Ys_ = bsxfun(@rdivide, Ys_, sig);
    tensorDenoiseWrapper(Ys_, r1, method, dataset);
  case 3
    load('/mat-files/DataMatched.mat');
    Ys_ = tensorDenoiseStandardize(DataMatched.Ys_);
    tensorDenoiseWrapper(Ys_, r1, method, dataset);
end

end
