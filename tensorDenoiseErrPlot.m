function [ summary ] = tensorDenoiseErrPlot( Y, trialCount )
%TENSORDENOISEERR2 function for mse vs err plots
% trial count is nxc matrix of trial counts


[n,t,c,r] = size(Y);
maxr = 10;

Ym = mean(Y,4,'omitnan');

errMeasure = 1;
perNeuron = 1;

% strategy 1 -- 
for rr = 1:maxr
  disp(num2str(rr));
  for ii = 1:10
    trial = [];
    Ys = zeros(n,t,c,rr);
    for nn = 1:n
    for cc = 1:c
      trial{nn,cc} = randi(trialCount(nn,cc),rr,1);
      Ys(nn,:,cc,:) = Y(nn,:,cc,trial{nn,cc});
    end
    end
  Ytr = mean(Ys,4,'omitnan');
  [mlRank,~] = tensorDenoiseRankEst( Ytr, 2);
  mse1(:,rr,ii) = tensorDenoiseERR(Ym, Ytr, errMeasure, perNeuron);
  for kk = 1:length(mlRank)
    lrModel = tensorDenoiseSVD(Ytr, mlRank{kk});
    mse2{kk}(:,rr,ii) = tensorDenoiseERR(Ym, lrModel, errMeasure, perNeuron);
  end
  end
end
summary.mse1 = squeeze(mse1);
summary.mse2 = squeeze(mse2);

end

 