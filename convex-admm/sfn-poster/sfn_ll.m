function ll = sfn_ll(Y,X)
% computes log likelihood
X = repmat(X,1,1,1,size(Y,4));
nn_inds = ~isnan(Y);

ll = -log(vec(poisspdf( Y(nn_inds), X(nn_inds) )));
ll = nanmean(ll);

% or

%ll = vec(poisspdf( Y(nn_inds), X(nn_inds) ));
%ll = nanmean(ll);
%log likelihood is good b/c we can take means instead of products...

end







