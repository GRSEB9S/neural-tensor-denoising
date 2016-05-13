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
  if ~isfield(options,'relative');
    options.relative = 0; end

  %% get options
  perNeuron = options.perNeuron;
  threshold = options.threshold;

  %% threshold data before computing error?
  if threshold
    Yhat = Yhat.*(Yhat>0);
  end
  
  %% divisive factor
  if options.relative && options.perNeuron
      Ydiv = sum(Y(:,:).^2,2);
  elseif options.relative
      Ydiv = norm(Y(:)).^2;
  else
      Ydiv = 1;
  end
  
  %% compute MSE
  Y = Y(:,:);
  Yhat = Yhat(:,:);
  if ~perNeuron    
    err = norm(Y(:)-Yhat(:))^2/Ydiv;
  else
    %err = norm(Y(:,:)-Yhat(:,:)).^2./Ydiv; % replace wtih sum(_^2,2)
    err = sum((Y(:,:)-Yhat(:,:)).^2,2)/Ydiv;
  end

  
  
  
end

