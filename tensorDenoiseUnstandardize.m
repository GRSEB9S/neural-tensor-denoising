function [ Y_out ] = tensorDenoiseUnstandardize( Y_in, mu, sigma2 )
%inverse of tensorDenoiseStandardize
% see section 2.6 of GLRM paper

  Y_out = bsxfun(@times, Y_in, sigma2);
  Y_out = bsxfun(@plus, Y_out, mu);
  
end

