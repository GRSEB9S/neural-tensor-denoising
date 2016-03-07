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

  %% threshold data?
  if threshold
    Yhat = Yhat.*(Yhat>0);
  end
  
  %% compute relative error
  Y = Y(:,:);
  Yhat = Yhat(:,:);
  if ~perNeuron    
    err = norm(Y(:)-Yhat(:))^2;%/norm(Y(:))^2);
  else
    err = norm(Y(:,:)-Yhat(:,:)).^2;%./norm(Y(:,:)).^2);
  end

  
  
  
end

