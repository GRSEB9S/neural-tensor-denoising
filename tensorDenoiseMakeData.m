function Data = tensorDenoiseMakeData()
%tensorDenoiseMakeData() - make data files for TD figures.

% If on HPC:
if isdir('/ifs/')
  cd /ifs/scratch/zmbbi/la_lab/jss2219/
  addpath(genpath('/ifs/scratch/zmbbi/la_lab/jss2219/'))
end

%% params
% standard deviation for Gaussian filter
% bin width for sampling
sd = 10;
bw = 10;

%% load Lara data
% Lara = preS('Datasets/Motor/Lara-20141105/Sstructs/', bw, sd);
% 
% Data = tensorizeDataStruct(Lara);
% Data.Ysm = mean(Data.Ys,4,'omitnan');
% Data.size = size(Data.Y);
% Data.trialCount = getTrialCount(Data.Y);

%% load Shenoy Data
Shenoy = preRC('RC,2009-09-18,1-2,good-ss', bw, sd);
aligned = 1;

Data.Y = Shenoy.aligned(aligned).Y;
Data.Ys = Shenoy.aligned(aligned).Ys;
Data.Y_long = Shenoy.aligned(aligned).Y_long;
Data.Ys_long = Shenoy.aligned(aligned).Ys_long;
Data.cond = Shenoy.cond;
clear Shenoy

Data = tensorizeDataStruct(Data);
Data.Ysm = mean(Data.Ys,4,'omitnan');
Data.size = size(Data.Y);
Data.trialCount = getTrialCount(Data.Y);

%% get sort indices
[~,trialSort] = sort(sum(Data.trialCount,2),'descend');
[~,rangeSort] = sort(range(Data.Ysm(:,:),2),'descend');
[~,varSort] = sort(var(Data.Ysm(:,:),[],2),'descend');
Ystd = std(Data.Ys,[],4,'omitnan'); % is this correct for snr?
[~,snrSort] = sort(range(Data.Ysm(:,:),2)./mean(Ystd(:,:),2)./sqrt(sum(Data.trialCount,2)),'descend');

%% sort neurons
sortType = varSort;

Data.Y = Data.Y(sortType,:,:,:);
Data.Ys = Data.Ys(sortType,:,:,:);
Data.Ysm = Data.Ysm(sortType,:,:);
Data.trialCount = Data.trialCount(sortType,:);
Data.sort = sortType;

%% select top neurons
maxn = 100;
Data.Y = Data.Y(1:maxn,:,:,:);
Data.Ys = Data.Ys(1:maxn,:,:,:);
Data.Ysm = Data.Ysm(1:maxn,:,:);

%% Normalize data (?) 
[Data.Ysm, div] = tensorDenoiseNormalize(Data.Ysm);
Data.Y = bsxfun(@rdivide, Data.Y, div);
Data.Ys = bsxfun(@rdivide, Data.Ys, div);
Data.div = div;

%% save
save('TensorDenoising/mat-files/Data.mat','Data');

end

