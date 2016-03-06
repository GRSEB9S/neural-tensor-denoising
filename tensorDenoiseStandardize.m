function [ Y_out, mu, sigma2 ] = tensorDenoiseStandardize( Y_in )
%Standardize neurons of tensor Y_in
% see section 2.6 of GLRM paper

  sigma2 = var(Y_in(:,:),[],2,'omitnan');
  mu = mean(Y_in(:,:),2,'omitnan');
  Y_out = bsxfun(@plus, Y_in, -mu);
  Y_out = bsxfun(@rdivide, Y_out, sigma2);

end

