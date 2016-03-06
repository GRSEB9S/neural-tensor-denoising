function [ rankEst, summary ] = tensorDenoiseRankEst( Ydata, method, values )
%TENSORDENOISERANKEST computes MSE for all possible ML ranks
% brute force -- speed up later.
% method 1 2 3
% values for method 2: thresholds
% values for method 3: multiplier

% e.g. make use of singular values?  
if nargin < 3
  values = [];
end

%%
N = ndims(Ydata);
size_tens = size(Ydata);

%%
switch method

  %% method 1 - smart grid search
  % maybe start with upper bounds
  % and then decrease lexicographically, moving onto the next mode once
  % hitting below some threshold. pick best based on
  % 1 number of elements in the model
  % 2 min core size, (r1+r2+r3)
  case 1
    
    
  %% method 2 - var explained thresholds on all modes
  case 2
    if isempty(values); ths = 0.95; end;
    for kk = 1:length(ths)
      [U,S,sv] = mlsvd(Ydata);
      loss = cellfun(@(i)cumsum(i(:).^2)/sum(i(:).^2), sv, 'UniformOutput', false);
      for nn = 1:N
        rankEst{kk}(nn) = find(loss{nn} > ths(kk),1);
      end
    end
    summary.thresholds = ths;
    
  %% method 3 - use mlrankest from tensorlab toolbox
  case 3
    if isempty(values); multiplier = logspace(-2,2,10); end;
    for kk = 1:length(multiplier)
      options.XMultiplier = multiplier(kk);
      rankEst{kk} = mlrankest( Ydata, options );
    end
    summary.multiplier = multiplier;
    
  %% method 4 - generalization?
    
end
   
% if length(rankEst) == 1
%   rankEst = rankEst{1};
% end

end

