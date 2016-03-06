function tensorDenoiseFiguresCluster(in)
% generate figure 2 data. for hpc cluster.

if nargin == 0
  in = 1:3;
end

%% load data from tensorDenoiseMakeData.m
load('Data.mat');
load('DataFull.mat');
load('DataMatched.mat');
load('DataSim.mat');

%% struct names
snames = {'real','sim','matched'};

%% options
options = [];
options.gridStep = 2;
options.minRank = [20 2 2];
options.maxRank = [80 20 20];
options.maxr = 20;
options.resample = 1;
options.verbose = 0;


%% real data
if ismember(1,in)
  Ys_ = tensorDenoiseStandardize(Data.Ys_);

  real(1).summary = tensorDenoiseGridSearchCV(Ys_, options );
  for mm = 1:3
    options.mode = mm; real(mm+1).summary = tensorDenoiseUnfoldingCV(Ys_, options );
  end
  %real(5).summary = tensorDenoiseFilterCV(Ys_, options );
end

%% simulated data
if ismember(2,in)
  [Ytrue, mu, sig] = tensorDenoiseStandardize(DataSim.Ygt);
  Ys_ = bsxfun(@plus, DataSim.Ys_, -mu);
  Ys_ = bsxfun(@rdivide, Ys_, sig);

  sim(1).summary = tensorDenoiseGridSearchCV( Ys_, options );
  for mm = 1:3
    options.mode = mm; sim(mm+1).summary = tensorDenoiseUnfoldingCV( Ys_, options );
  end
  %sim(5).summary = tensorDenoiseFilterCV( Ys_, options );
end
  
%% matched data
if ismember(3,in)
  Ys_ = tensorDenoiseStandardize(DataMatched.Ys_);

  matched(1).summary = tensorDenoiseGridSearchCV( Ys_, options );
  for mm = 1:3
    options.mode = mm; matched(mm+1).summary = tensorDenoiseUnfoldingCV( Ys_, options );
  end
  %matched(5).summary = tensorDenoiseFilterCV( Ys_, options );
end

%% save
save(['/ifs/home/zmbbi/la_lab/jss2219/clusterOutput-' datestr(now,'ddmmyy-HHMM')])

end
