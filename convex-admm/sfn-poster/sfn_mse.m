function mse = sfn_mse(Y,X)

   X = X.*(X>0);
   mse1 = vec(nanmean(Y,4) - X).^2;
   mse1 = nanmean(mse1);
   
   % or
   %X = repmat(X,1,1,1,size(Y,4));
   %nn_inds = ~isnan(Y);
   %mse2 = vec(Y(nn_inds) - X(nn_inds)).^2;
   %mse2 = nanmean(mse2);
   
   mse = mse1;
   

end



