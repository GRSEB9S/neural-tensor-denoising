function [ Y_out, mu, sigma2 ] = tensorDenoiseStandardize( Y_in, mu, sigma2 )
%Standardize neurons of tensor Y_in
% see section 2.6 of GLRM paper

  if nargin == 1
    sigma2 = var(Y_in(:,:),[],2,'omitnan');
    mu = mean(Y_in(:,:),2,'omitnan');
  end
  
  Y_out = bsxfun(@plus, Y_in, -mu);
  Y_out = bsxfun(@rdivide, Y_out, sigma2);

end

