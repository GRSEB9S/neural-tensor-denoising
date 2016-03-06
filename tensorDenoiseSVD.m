function [ Ysvd, size_core ] = tensorDenoiseSVD( Y, threshold )
%TENSORDENOISESVD truncated multilinear svd on Y
% threshold - scalar in [0, 1] for estimating size_core
% threshold - vector for manually setting size_core

  if isscalar(threshold)
    [~,~,sv] = mlsvd(Y);
    loss = cellfun(@(i)cumsum(i(:).^2)/sum(i(:).^2), sv, 'UniformOutput', false);
    for nn = 1:ndims(Y)
      size_core(nn) = find(loss{nn} > threshold, 1);
    end
  else
    size_core = threshold;
  end
  
  [U,S] = mlsvd(Y, size_core);
  Ysvd = lmlragen(U,S);
   
end

