function [ summary ] = tensorDenoiseGridSearchCV( Y, options )
% grid search through core tensor sizes.
% options.ths - SVD threshold for setting minRank and maxRank
% options.

  %% set defaults for options structure
  if nargin < 2
    options = struct; end
  if ~isfield(options,'coreMeasure');
    options.coreMeasure = 2; end
  if ~isfield(options,'gridStep');
    options.gridStep = 1; end
  if ~isfield(options,'r1');
    options.r1 = 10; end
  if ~isfield(options, 'resample'); % run resampleTrials()?
    options.resample = 0; end
  if ~isfield(options, 'verbose'); 
    options.verbose = 1; end
  if ~isfield(options, 'rng_iter'); 
    options.rng_iter = 0; end
  if ~isfield(options, 'method'); 
    options.method = 0; end

  %% get options
  gridStep = options.gridStep;
  r1 = options.r1;
  resample = options.resample;
  verbose = options.verbose;
  rng_iter = options.rng_iter;
  method = options.method;
  minRank = options.minRank;
  maxRank = options.maxRank;
  
  %% define set of multilinear ranks to search over
  Ym = mean(Y,4,'omitnan'); 
  N = ndims(Ym);
  ysize = size(Ym);
  
  if method == 0
    in = arrayfun(@(i)minRank(i):gridStep:maxRank(i),1:N,'UniformOutput',false);
    out = cell(1,N);
    [out{:}] = ndgrid(in{:});
    out = cell2mat(cellfun(@(i)i(:),out,'UniformOutput',false)); 
  else
    out = repmat(ysize, ysize(method), 1);
    out(:,method) = 1:ysize(method);
  end
  
  %% core complexity measures
  model_elts = @(s)(sum(bsxfun(@times,size(Ym),s),2)+prod(s,2));
  core_elts = @(s)(prod(s,2));
  core_sum = @(s)(sum(s,2));

  %% resample trials
  if resample
    Y = resampleTrials(Y, 1, r1+1000*rng_iter);
  end
  
  %% compute error for each multilinear rank
  erroptions = struct;
  erroptions.threshold = 1;
  erroptions.perNeuron = 0;
  err = zeros(size(out,1), r1);

  if verbose; disp('r2 / r1'); end
  
  for r2 = 1:r1
    if verbose; disp([num2str(r2) ' / ' num2str(r1)]); end
    Ytrain = Y(:,:,:,1:r1);
    Ytrain(:,:,:,r2) = [];
    Yval = Y(:,:,:,r2);
    for ii = 1:size(out,1)     
      % compute error:
      Yhat = tensorDenoiseSVD(mean(Ytrain, 4, 'omitnan'), out(ii,:));
      err(ii,r2) = tensorDenoiseERR( Yval, Yhat, erroptions );
    end
  end
  err_avg = mean(err,2);

  
  %% optimal model -- lots to do
  [p_out, p_complexity, p_err_avg] = paretoPoints(out, err_avg, core_elts);
  
  % search through possible slope values:
  % f(x) - mx. find min for different choices of m.
  
  minind = find(err_avg == min(err_avg),1);
  minrank = out(minind,:);
  
  mrange = 0:100:10000;
  for mm = 1:length(mrange)
    L = err_avg + mrange(mm)*core_elts(out);
    minind_p(mm,1) = find(L == min(L),1);
    minrank_p(mm,:) = out(minind_p(mm),:);
  end
  
  %Yest = tensorDenoiseSVD(Ysm, minrank);
  
  %% output 
  summary.ranks = out;
  summary.model_elts = model_elts(out);
  summary.core_elts = core_elts(out);
  summary.core_sum = core_sum(out);
  summary.err = err_avg;
  summary.p_err = p_err_avg;
  summary.p_complexity = p_complexity;
  summary.p_ranks = p_out;
  summary.minind = minind;
  summary.minrank = minrank;
  summary.minind_p = minind_p;
  summary.minrank_p = minrank_p;
  summary.options = options;
  

end