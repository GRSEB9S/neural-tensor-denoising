function DataSim = tensorDenoiseMakeDataSim(Data, size_core)

  %% basic params
  bw = 10;
  sd = 10;
  [n,t,c,r] = size(Data.Ys);
  
  %% get ground truth data
  Ygt = tensorDenoiseSVD(Data.Ysm, size_core);
  %Ygt = filterGauss(Ygt, 2, 2);
  Ygt = Ygt.*(Ygt>0);

  
  % expand in time
  YgtFull = permute(Ygt,[2 1 3]);
  YgtFull = 1/bw*interp1(1:size(YgtFull,1), YgtFull, linspace(1,size(YgtFull,1),bw*size(YgtFull,1)), 'linear', 'extrap');
  YgtFull = permute(YgtFull,[2 1 3]);
  YgtFull = YgtFull.*(YgtFull>0);

  % simulate spikes
  r_sim = 20;
  rng(42);
  DataSim.Y = poissrnd(mean(1/bw*Data.Y(:),'omitnan')/mean(YgtFull(:))*repmat(YgtFull,1,1,1,r_sim));
  DataSim.Ys = filterGauss(DataSim.Y, sd, 2);
  DataSim.Ysm = mean(DataSim.Ys,4,'omitnan');

  % spike binning
  DataSim.Y = permute(DataSim.Y, [1 3 4 2]);
  DataSim.Y = sum( reshape(DataSim.Y, [n c r_sim t bw]), 5);
  DataSim.Y = permute(DataSim.Y, [1 4 2 3]);

  % subsample
  DataSim.Ys = bw*DataSim.Ys(:,1+round(0.5*bw):bw:end,:,:);
  DataSim.Ysm = bw*DataSim.Ysm(:,1+round(0.5*bw):bw:end,:);
  DataSim.Ygt = bw*YgtFull(:,1+round(0.5*bw):bw:end,:);
  DataSim.trialCount = r_sim*ones(n,c);
  DataSim.r_sim = r_sim;

end