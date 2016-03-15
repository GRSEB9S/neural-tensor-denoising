function [ err, err2 ] = tensorDenoiseERR( Y, Yhat, options )
%TENSORDENOISEMSE Computes MSE of tensors Y and X
%   Y is data, Yhat is estimate of data.

  %% options struct
  if nargin < 3
    options = struct; end
  
  if ~isfield(options,'perNeuron');
    options.perNeuron = 0; end
  if ~isfield(options,'errMeasure');
    options.errMeasure = 1; end
  if ~isfield(options,'threshold');
    options.threshold = 1; end

  %% get options
  perNeuron = options.perNeuron;
  threshold = options.threshold;

  %% threshold data before computing error?
  if threshold
    Yhat = Yhat.*(Yhat>0);
  end
  
  %% compute relative error
  Y = Y(:,:);
  Yhat = Yhat(:,:);
  if ~perNeuron    
    err = 1/length(Y(:))*norm(Y(:)-Yhat(:))^2;
  else
    err = 1/size(Y(:,:),2)*norm(Y(:,:)-Yhat(:,:)).^2;
  end

  
  
  
end

