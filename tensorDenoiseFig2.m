function outStruct = tensorDenoiseFig2(Data)

% implement: bootstrap
options = [];
options.gridStep = 3;
options.minRank = [20 6 6];
options.maxRank = [60 20 20];
options.maxr = 10;
options.resample = 1;

outStruct(1).summary = tensorDenoiseGridSearchCV( tensorDenoiseStandardize(Data.Ys_), options );
for mm = 1:3
  options.mode = mm; outStruct(mm+1).summary = tensorDenoiseUnfoldingCV( tensorDenoiseStandardize(Data.Ys_), options );
end
outStruct(5).summary = tensorDenoiseFilterCV( tensorDenoiseStandardize(Data.Ys_), options );


%% plot results
figure; hold all

options = struct;
options.threshold = 0;

Ytrue = mean(tensorDenoiseStandardize(Data.Ys_),4,'omitnan');
for rr = 1:length(outStruct(1).summary.errs)
  Ytr = mean(tensorDenoiseStandardize(Data.Ys_(:,:,:,1:rr)), 4, 'omitnan');
  errTrue = tensorDenoiseERR( Ytrue, Ytr, options);
  errDN  = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, outStruct(1).summary.minrank(rr,:)), options);
  errDN1 = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, outStruct(2).summary.minrank(rr,:)), options);
  errDN2 = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, outStruct(3).summary.minrank(rr,:)), options);
  errDN3 = tensorDenoiseERR( Ytrue, tensorDenoiseSVD(Ytr, outStruct(4).summary.minrank(rr,:)), options);
  errFilt = tensorDenoiseERR( Ytrue, filterGauss(Ytr, outStruct(5).summary.minfilt(rr), 2), options); 
  
  vaf(rr,1) = errTrue.vaf;
  vaf(rr,2) = errDN.vaf;
  vaf(rr,3) = errDN1.vaf;
  vaf(rr,4) = errDN2.vaf;
  vaf(rr,5) = errDN3.vaf;
  vaf(rr,6) = errFilt.vaf;
end

figure; hold all
plot(vaf)
xlabel('trials (per neuron per condition)');
ylabel('1 - VAF');

end