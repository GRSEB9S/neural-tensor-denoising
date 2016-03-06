function tensorDenoiseSummaryPlot(summary_in, Ytrue, Ys_, size_core)

if nargin < 4
  size_core = [];
end

options = struct;
options.threshold = 0;

for rr = 1:length(summary_in(1).summary.errs)
  Ytr = mean(Ys_(:,:,:,1:rr), 4, 'omitnan');
  errTrue = tensorDenoiseERR( Ytrue, Ytr, options);
  errDN  = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, summary_in(1).summary.minrank(rr,:)), options);
  errDN1 = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, summary_in(2).summary.minrank(rr,:)), options);
  errDN2 = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, summary_in(3).summary.minrank(rr,:)), options);
  errDN3 = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, summary_in(4).summary.minrank(rr,:)), options);
  errFilt = tensorDenoiseERR( Ytrue, filterGauss(Ytr, summary_in(5).summary.minfilt(rr), 2), options); 
  
  vaf(rr,1) = errTrue.vaf;
  vaf(rr,2) = errDN.vaf;
  vaf(rr,3) = errDN1.vaf;
  vaf(rr,4) = errDN2.vaf;
  vaf(rr,5) = errDN3.vaf;
  vaf(rr,6) = errFilt.vaf;
end

figure; hold all
xaxis = 2:length(summary_in(1).summary.errs)+1;
plot(xaxis,vaf)
legend('average','tensor','neuron','time','condition','filter');
xlabel('trials (per neuron per condition)');
ylabel('1 - VAF, relative to "ground truth" data');
title('main analysis | relative error plot')

figure; hold all
plot(xaxis,summary_in(1).summary.minrank);
legend('neuron','time','condition');
if ~isempty(size_core)
  plot([xaxis(1) xaxis(end)],[size_core; size_core],'--');
  legend('neuron','time','condition','true neuron','true time','true condition');
end
xlabel('trials (per neuron per condition)');
ylabel('rank');
title('optimal multilinear rank | tensor denoising');

figure; hold all
plot(xaxis,summary_in(2).summary.minrank(:,1));
plot(xaxis,summary_in(3).summary.minrank(:,2));
plot(xaxis,summary_in(4).summary.minrank(:,3));
legend('neuron','time','condition');
if ~isempty(size_core)
  plot([xaxis(1) xaxis(end)],[size_core; size_core],'--');
  legend('neuron','time','condition','true neuron','true time','true condition');
end
xlabel('trials (per neuron per condition)');
ylabel('rank');
title('optimal ranks | matrix denoising')

figure; hold all
plot(xaxis,summary_in(5).summary.minfilt);
legend('filter width');
xlabel('trials (per neuron per condition)');
ylabel('filter width (SD for Gaussian filter)');
title('optimal filter width | filter denoising');

end
