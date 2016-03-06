function [ ranks, complexity, err ] = tensorDenoiseGridSearch( Y, options )
% grid search through core tensor sizes.

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
  
  %% get options
  ths = options.ths;
  coreMeasure = options.coreMeasure;
  errMeasure = options.errMeasure;
  gridStep = options.gridStep;
  
  %% svd / define grid search box
  N = ndims(Y);
  [U,S,sv] = mlsvd(Y);
  loss = cellfun(@(i)cumsum(i(:).^2)/sum(i(:).^2), sv, 'UniformOutput', false);
  minRank = zeros(N,1);
  maxRank = zeros(N,1);
  for nn = 1:N
    minRank(nn) = find(loss{nn} > ths(1), 1);
    maxRank(nn) = find(loss{nn} > ths(2), 1);
  end
  
  %% grid
  in = arrayfun(@(i)minRank(i):gridStep:maxRank(i),1:N,'UniformOutput',false);
  out = cell(1,N);
  [out{:}] = ndgrid(in{:});
  out = cell2mat(cellfun(@(i)i(:),out,'UniformOutput',false));
  
  %% error
  err = zeros(size(out,1), 1);
  for ii = 1:size(out,1)
    UU = cell(N,1);
    sinds = cell(N,1);
    for nn = 1:N
      UU{nn} = U{nn}(:,1:out(ii,nn));
      sinds{nn} = 1:out(ii,nn);
    end
    SS = S(sinds{:});
    % compute error:
    err(ii) = tensorDenoiseERR( Y, lmlragen(UU, SS), errMeasure );
    %disp([num2str(ii) ' / ' num2str(size(out,1))]);
  end

  %% core complexity measures
  model_elts = @(s)(sum(bsxfun(@times,size(Y),s),2)+prod(s,2));
  core_elts = @(s)(prod(s,2));
  core_sum = @(s)(sum(s,2));
  
  switch coreMeasure
    case 1; mem = model_elts;
    case 2; mem = core_elts;
    case 3; mem = core_sum;
  end
  
  %% pareto optimal
  x = mem(out);
  y = err;
  
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
  
  %% output
  ranks = out;
  complexity = x;
  err = y;
    
end
