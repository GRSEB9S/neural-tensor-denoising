function [ Y_out ] = tensorDenoiseNormalize( Y_in, const )
% normalize non-neagative neurons of tensor Y_in

  if nargin < 2
    const = 5;
  end
  Y_out = bsxfun(@rdivide, Y_in, const + range(Y_in(:,:),2));
end

