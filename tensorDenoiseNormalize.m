function [ Y_out, div ] = tensorDenoiseNormalize( Y_in, const )
% normalize non-neagative neurons of tensor Y_in

  if nargin < 2
    const = 0.5*max(Y_in(:));
  end
  div = max(Y_in(:,:),[],2) + const;
  Y_out = bsxfun(@rdivide, Y_in, div);
  
end

