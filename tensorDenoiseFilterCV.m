function [ summary ] = tensorDenoiseFilterCV( Y, options )
% sweep across filter widths and cross validate.
% gaussian smoothing, etc...
% NOT FINISHED



  %% set defaults for options structure
  if nargin < 2
    options = struct; end
  if ~isfield(options,'maxr');
    options.maxr = 10; end
  if ~isfield(options, 'resampleTrials'); % for mixing up trial counts
    options.resample = 0; end
  if ~isfield(options, 'mode');
    options.mode = 1; end
  
  %% LOOCV options  
  [n,t,c,r] = size(Y);
  Ys = Y;
  Ysm = mean(Ys,4,'omitnan');
  
  %% get options
  maxr = options.maxr;
  resample = options.resample;

  %% grid
  out = 0.5:0.25:7;

  %% error
  options = struct;
  %% edit
  options.threshold = 0; % usually set to 1
  %% /edit
  options.perNeuron = 0;
  
  disp(['r2 / r1 / maxr']);

  for r1 = 2:maxr
    % replacement?
    if resample
      rng(r1);
      Ys = resampleTrials(Y, 1);
    end
    
    for r2 = 1:r1
      disp([num2str(r2) ' / ' num2str(r1) ' / ' num2str(maxr)]);
      Ytrain = Ys(:,:,:,1:r1);
      Ytrain(:,:,:,r2) = [];
      Yval = Ys(:,:,:,r2);
      for ii = 1:length(out)     
        % compute error:
        err = tensorDenoiseERR( Yval, mean(filterGauss(Ytrain,out(ii),2),4,'omitnan'), options );        
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
    minfilt(r1-1) = out(minind);
    
  end
  
  %% output
  summary.filts = out;
  summary.errs = err_summary;
  summary.minfilt = minfilt;
  
  summary.options = options;
  
  %% function not ready
  function [ranks,x,y] = pareto(out, coreMeasure)
    %% pareto optimal
    out_old = out;
    for r1 = 2:maxr
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