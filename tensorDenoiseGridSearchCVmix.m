function [ summary ] = tensorDenoiseGridSearchCVmix( Y, options )
% runs tensorDenoiseGridSearchCV but over several instantiations of mixing
% / replacing trials


  %% set defaults for options structure
  % options: with replacement or not
  
  if nargin < 2
    options = struct; end
  if ~isfield(options, 'replacement')
    options.replacement = 1; end
  if ~isfield(options, 'numRuns')
    options.numRuns = 10; end
  
  %% mix up trials
  cvoptions = [];
  cvoptions.ths = [0.7 0.999];
  cvoptions.coreMeasure = 2;
  cvoptions.errMeasure = 1;
  cvoptions.gridStep = 5;
  cvoptions.minRank = [8 3 3];
  cvoptions.maxRank = [80 30 24];
  for rr = 1:options.numRuns
   disp(['bootstrap # ' num2str(rr)]);
    summary{rr} = tensorDenoiseGridSearchCV( Y, cvoptions );
  end
  
end






