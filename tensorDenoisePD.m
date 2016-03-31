%% compute preferred direction
% take in (n,t,c) tensor.
% for each (n,t), compute the PD.
% possibly consider longer time windows, e.g. sliding time bin averaging.

Y = Data.Ysm;
[n,t,c] = size(Y);

%% Load target data
load('/Users/Jeff/Documents/MATLAB/Datasets/Motor/Lara-20141105/NB.mat');

%% get targ values;
% note: each N(1), N(2), etc, has a different x,y, coordinate for the
% targets, but seem to be essentially at the same angle. so as long as we
% do arctan we should be fine(?)
targ = zeros(2,c);
for cc = 1:c
  targ(1,cc) = N(1).condTable(cc,5).TARGINFO.x/150; % put in sensible range /150.
  targ(2,cc) = N(1).condTable(cc,5).TARGINFO.y/150;
end

%% get PD
pd = zeros(n,t);
for nn = 1:n
  for tt = 1:t
    tc = squeeze(Y(nn,tt,:))';
    pvec = tc/targ;
    pd(nn,tt) = atan2(pvec(1),pvec(2))*180/pi;
  end
end

%% get sliding window PD
Ys = filterGauss(Y,5,2);
pd2 = zeros(n,t);
for nn = 1:n
  for tt = 1:t
    tc = squeeze(Ys(nn,tt,:))';
    pvec = tc/targ;
    pd2(nn,tt) = atan2(pvec(1),pvec(2))*180/pi;
  end
end

%% plots
figure; hold all;
pdplot = cat(3,pd,pd2);
tensorSubplot(pdplot);










