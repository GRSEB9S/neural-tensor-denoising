function [ Data, DataFull, DataSim, DataMatched ] = tensorDenoiseMakeData()
%tensorDenoiseMakeData() - make data files for TD figures.

cd /ifs/scratch/zmbbi/la_lab/jss2219/
addpath(genpath('/ifs/scratch/zmbbi/la_lab/jss2219/'))

%% params
% standard deviation for Gaussian filter
% bin width for sampling
sd = 10;
bw = 10;

%% load data
% Lara = preS2('Datasets/Motor/Lara-20141105/Sstructs/', bw, sd);
Shenoy = preRC2('RC,2009-09-18,1-2,good-ss', bw, sd);
Shenoy.Y = Shenoy.aligned(4).A;
Shenoy.Ys = Shenoy.aligned(4).As;

Data = compressTrials2(Shenoy);
Data.Ysm = mean(Data.Ys,4,'omitnan');
Data.size = size(Data.Y);
Data.trialCount = getTrialCount(Data.Y);
% for nn = 1:Data.size(1)
%   Data.trialCount(nn,:) = histc(Data.cond{nn},1:Data.size(3));
% end

[n t c r] = size(Data.Y);

% *Full = 1 time bin = 1 ms
% LaraFull = preS2('Datasets/Motor/Lara-20141105/Sstructs/', 1, sd);
ShenoyFull = preRC2('mat-files/RC,2009-09-18,1-2,good-ss', 1, sd);
ShenoyFull.Y = Shenoy.aligned(4).A;
ShenoyFull.Ys = Shenoy.aligned(4).As;

DataFull = compressTrials2(ShenoyFull);
DataFull.Ysm = mean(DataFull.Ys, 4, 'omitnan');
DataFull.size = size(DataFull.Y);
Data.trialCount = getTrialCount(DataFull.Y);
% for nn = 1:DataFull.size(1)
%   DataFull.trialCount(nn,:) = histc(DataFull.cond{nn},1:DataFull.size(3));
% end

%% get sort indices
% trial counts. maybe just implement into compresstrials2
% take best neurons
[~,trialSort] = sort(sum(Data.trialCount,2),'descend');
[~,rangeSort] = sort(range(Data.Ysm(:,:),2),'descend');
Ystd = std(Data.Ys,[],4,'omitnan'); % is this correct for snr?
[~,varSort] = sort(var(Data.Ys(:,:)',1),'descend');
[~,snrSort] = sort(range(Data.Ysm(:,:),2)./mean(Ystd(:,:),2)./sqrt(sum(Data.trialCount,2)),'descend');

%% sort neurons
sortType = varSort;

Data.Y = Data.Y(sortType,:,:,:);
Data.Ys = Data.Ys(sortType,:,:,:);
Data.Ysm = Data.Ysm(sortType,:,:);
Data.trialCount = Data.trialCount(sortType,:);
Data.sort = sortType;

DataFull.Y = DataFull.Y(sortType,:,:,:);
DataFull.Ys = DataFull.Ys(sortType,:,:,:);
DataFull.Ysm = DataFull.Ysm(sortType,:,:);
DataFull.trialCount = DataFull.trialCount(sortType,:);
DataFull.sort = sortType;

%% select top 60 only (or top 80?)
maxn = 100;
Data.Y_ = Data.Y(1:maxn,:,:,:);
Data.Ys_ = Data.Ys(1:maxn,:,:,:);
Data.Ysm_ = Data.Ysm(1:maxn,:,:);

DataFull.Y_ = DataFull.Y(1:maxn,:,:,:);
DataFull.Ys_ = DataFull.Ys(1:maxn,:,:,:);
DataFull.Ysm_ = DataFull.Ysm(1:maxn,:,:,:);

%% simulated data
DataSim.th = 0.90;
% get size_core from standardized data
[Y_standard, mu, sig] = tensorDenoiseStandardize(Data.Ysm_);
[Ygt,DataSim.size_core] = tensorDenoiseSVD(Y_standard, DataSim.th);
Ygt = tensorDenoiseUnstandardize(Ygt, mu, sig);
Ygt = filterGauss(Ygt, 2, 2);
DataSim.Ygt = Ygt.*(Ygt>0);

% interpolate data
YgtFull = permute(DataSim.Ygt,[2 1 3]);
YgtFull = interp1(1:size(YgtFull,1), YgtFull, linspace(1,size(YgtFull,1),bw*size(YgtFull,1)), 'linear', 'extrap');
YgtFull = permute(YgtFull,[2 1 3]);
YgtFull = YgtFull.*(YgtFull>0);

r_sim = 40;
rng(42)
DataSim.Y_ = poissrnd(mean(DataFull.Y(:),'omitnan')/mean(YgtFull(:))*repmat(YgtFull,1,1,1,r_sim));
DataSim.Ys_ = filterGauss(DataSim.Y_, sd, 2);
DataSim.Ysm_ = mean(DataSim.Ys_,4,'omitnan');

% write function to do spike binning this way.
DataSim.Y_ = permute(DataSim.Y_, [1 3 4 2]);
DataSim.Y_ = sum( reshape(DataSim.Y_, [maxn c r_sim t bw]), 5);
DataSim.Y_ = permute(DataSim.Y_, [1 4 2 3]);

% DataSim.Ygt = DataSim.Ygt(:,1:bw:end,:);
DataSim.Ys_ = DataSim.Ys_(:,1:bw:end,:,:);
DataSim.Ysm_ = DataSim.Ysm_(:,1:bw:end,:);
DataSim.trialCount = r_sim*ones(n,c);
DataSim.r_sim = r_sim;

%% simulated data - 'fake-real'
DataMatched = DataSim;
DataMatched.Y_ = nan([maxn t c r]);
DataMatched.Ys_ = nan([maxn t c r]);

for nn = 1:maxn
  for cc = 1:c   
    DataMatched.Y_(:,:,:,1:Data.trialCount(nn,cc))   = DataSim.Y_(:,:,:,1:Data.trialCount(nn,cc));
    DataMatched.Ys_(:,:,:,1:Data.trialCount(nn,cc))  = DataSim.Ys_(:,:,:,1:Data.trialCount(nn,cc));
  end
end

DataMatched.Ysm_ = mean(DataMatched.Ys_,4,'omitnan');

%% save
save('mat-files/Data_maze.mat','Data');
save('mat-files/DataFull_maze.mat','DataFull');
save('mat-files/DataSim_maze.mat','DataSim');
save('mat-files/DataMatched_maze.mat','DataMatched');

end

