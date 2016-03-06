function [ summary ] = tensorDenoiseGridSearchCV( Y, options )
% grid search through core tensor sizes.
% options.ths - SVD threshold for setting minRank and maxRank
% options.

  %% set defaults for options structure
  if nargin < 2
    options = struct; end
  if ~isfield(options,'ths');
    options.ths = [ 0.95 0.99 ]; end
  if ~isfield(options,'coreMeasure');
    options.coreMeasure = 2; end
  if ~isfield(options,'errMeasure');
    options.errMeasure = 1; end
  if ~isfield(options,'gridStep');
    options.gridStep = 1; end
  if ~isfield(options,'r_range');
    options.r_range = 2:10; end
  if ~isfield(options, 'resampleTrials'); % for mixing up trial counts
    options.resample = 0; end
  if ~isfield(options, 'verbose');
    options.verbose = 1; end

  %% LOOCV options          
  [n,t,c,r] = size(Y);
  Ys = Y;
  Ysm = mean(Ys,4,'omitnan');
  
  %% get options
  ths = options.ths;
  gridStep = options.gridStep;
  r_range = options.r_range;
  resample = options.resample;
  verbose = options.verbose;
  
  %% svd / define grid search box
  N = ndims(Ysm);
  [U,S,sv] = mlsvd(Ysm);
  loss = cellfun(@(i)cumsum(i(:).^2)/sum(i(:).^2), sv, 'UniformOutput', false);
  minRank = zeros(N,1);
  maxRank = zeros(N,1);
  for nn = 1:N
    minRank(nn) = find(loss{nn} > ths(1), 1);
    maxRank(nn) = find(loss{nn} > ths(2), 1);
  end
  % or overwrite if specified
  if isfield(options,'minRank');
    minRank = options.minRank; end
  if isfield(options,'maxRank');
    maxRank = options.maxRank; end
  
  %% grid
  in = arrayfun(@(i)minRank(i):gridStep:maxRank(i),1:N,'UniformOutput',false);
  out = cell(1,N);
  [out{:}] = ndgrid(in{:});
  out = cell2mat(cellfun(@(i)i(:),out,'UniformOutput',false));
  
  %% core complexity measures
  model_elts = @(s)(sum(bsxfun(@times,size(Ysm),s),2)+prod(s,2));
  core_elts = @(s)(prod(s,2));
  core_sum = @(s)(sum(s,2));

  %% error
  options = struct;
  %% edit
  options.threshold = 0; % usually set to 1
  %% /edit
  options.perNeuron = 0;
  
  if verbose;
  disp(['r2 / r1 / maxr']); end

  for r1 = r_range
    % replacement?
    if resample
      rng(r1);
      Ys = resampleTrials(Y, 1);
    end
    
    for r2 = 1:r1
      if verbose
      disp([num2str(r2) ' / ' num2str(r1) ' / ' num2str(r_range(end))]); end
      Ytrain = Ys(:,:,:,1:r1);
      Ytrain(:,:,:,r2) = [];
      Yval = Ys(:,:,:,r2);
      for ii = 1:size(out,1)     
        % compute error:
        err = tensorDenoiseERR( Yval, tensorDenoiseSVD(mean(Ytrain, 4, 'omitnan'), out(ii,:)), options );        
        err_summary(r1-1).mse(ii,r2) = err.mse;
        err_summary(r1-1).rel(ii,r2) = err.rel;
        err_summary(r1-1).vaf(ii,r2) = err.vaf;
      end
    end
    err_summary(r1-1).mse_avg = mean(err_summary(r1-1).mse,2);
    err_summary(r1-1).rel_avg = mean(err_summary(r1-1).rel,2);
    err_summary(r1-1).vaf_avg = mean(err_summary(r1-1).vaf,2);
    
    % get mins
    minind = find(err_summary(r1-1).mse_avg == min(err_summary(r1-1).mse_avg),1);
    minrank(r1-1,:) = out(minind,:);

  end
  
  %% output
  summary.ranks = out;    
  summary.model_elts = model_elts(out);
  summary.core_elts = core_elts(out);
  summary.core_sum = core_sum(out);
  summary.errs = err_summary;  
  summary.minrank = minrank;
  
  summary.options = options;
  


  %% function not ready
  function [ranks,x,y] = pareto(out, coreMeasure)
    %% pareto optimal
    out_old = out;
    for r1 = 2:r_range
      out = out_old;
      x = mem(out);
      y = err_summary.mse_avg(:,r1-1);

      d = sqrt(x.^2+y.^2);
      [d,idx] = sort(d);
      out = out(idx,:);
      x = x(idx);
      y = y(idx);

      p = 1;
      while sum(p) <= length(x) % originally it was ~=. now it is <=
          idx = (x >= x(p) & y > y(p)) | (x > x(p) & y >= y(p));
          if any(idx)
              out = out(~idx,:); x = x(~idx); y = y(~idx); d = d(~idx);
          end
          p = p+1-sum(idx(1:p));
      end
      [~,idx] = sortrows([mem(out) -y]);
      out = out(idx,:); x = x(idx); y = y(idx);
    end
  end

  %% resample trials
  % to do
  function Yout = resampleTrials( Yin, repl )
      trialCount = getTrialCount(Yin);
      for nn = 1:n;
        for cc = 1:c
          % replacement option 1: with replacement
          if repl
            trialInds = randi(trialCount(nn,cc), 1, trialCount(nn,cc));
          % replacement option 2: without replacement
          else
            trialInds = randperm(trialCount(nn,cc));
          end
          Yin(nn,:,cc,1:trialCount(nn,cc)) = Yin(nn,:,cc,trialInds);
        end
      end
      Yout = Yin;
  end

  %% get trial count. used in resampleTrials()
  function trialOut = getTrialCount(Yin)
    [n,t,c,r] = size(Yin);
    for nn = 1:n
      for cc = 1:c
        trialOut(nn,cc) = sum(~isnan(Yin(nn,1,cc,:)));
      end
    end
  end


end