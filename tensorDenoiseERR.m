function [ err ] = tensorDenoiseERR( Y, Yhat, options )
%TENSORDENOISEMSE Computes Error of tensors Y and X
%   Y is data, Yhat is estimate of data.

  %% options struct
  if nargin < 3
    options = struct; end
  
  if ~isfield(options,'perNeuron');
    options.perNeuron = 0; end
  if ~isfield(options,'errMeasure');
    options.errMeasure = 1; end
  if ~isfield(options,'threshold');
    options.threshold = 0; end

  %% get options
  perNeuron = options.perNeuron;
  threshold = options.threshold;

  %% threshold data before computing error?
  if threshold
    Yhat = Yhat.*(Yhat>0);
  end
  
  %% compute MSE
  Y = Y(:,:);
  Yhat = Yhat(:,:);
  if ~perNeuron    
    err = norm(Y(:)-Yhat(:))^2;
  else
    err = norm(Y(:,:)-Yhat(:,:)).^2; % replace wtih sum(_^2,2)
  end

  
  
  
end

