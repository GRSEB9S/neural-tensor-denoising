function [ err, err2 ] = tensorDenoiseERR( Y, Yhat, options )
%TENSORDENOISEMSE Computes MSE of tensors Y and X
%   Y is data, Yhat is estimate of data.
%   standard to use Y as spiking data or filtered data?
%   in sfn poster it was spiking data.

% measure 1 - MSE
% measure 2 - rel error: frob(Y-Yhat)/frob(Y)
% measure 3 - variance unexplained (1-VAF)
% measure 4 - lmlraerr.m? subspace angles

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
  
%% old approach
%   %% 
%   if ~perNeuron
%     switch errMeasure
%       case 1 % MSE
%         err = (Y - Yhat).^2;
%         err = mean(err(:),'omitnan');
%       case 2 % Rel Error
%         err = norm(Y(:) - Yhat(:),'fro')/norm(Y(:),'fro'); % cant handle nans
%       case 3 % 1 - VAF
%         err = min(1,norm(Y(:)-Yhat(:)).^2/norm(Y(:)).^2); % cant handle nans
%       case 4 % 
%     end
%   else
%     switch errMeasure
%       case 1
%         err = (Y - Yhat).^2;
%         err = mean(err(:,:),2,'omitnan');
%       case 2
%         for nn = 1:size(Y,1)      
%           err(nn,1) = norm(Y(nn,:) - Yhat(nn,:),'fro')/norm(Y(nn,:),'fro');
%         end
%       case 3
%         err = min(1,norm(Y(:,:)-Yhat(:,:)).^2./norm(Y(:,:)).^2); % double check to see if it works
%     end
%   end

  %% alternative approach
  err = struct;
  
  Y = Y(:,:);
  Yhat = Yhat(:,:);
  
  if ~perNeuron
    
    err.mse = mean((Y(:) - Yhat(:)).^2,'omitnan');
    err.rel = norm(Y(:) - Yhat(:),'fro')/norm(Y(:),'fro');
    err.vaf = min(1, norm(Y(:) - Yhat(:)).^2/norm(Y(:)).^2);
  
  else
    
    err2.mse = mean((Y(:,:)-Yhat(:,:)).^2,2,'omitnan');
    for nn = 1:size(Y,1)
      err2.rel(nn,1) = norm(Y(nn,:) - Yhat(nn,:), 'fro')/norm(Y(nn,:),'fro'); % needs fixing. nans -> []
    end
    err2.vaf(nn,1) = min(1,norm(Y(:,:)-Yhat(:,:)).^2./norm(Y(:,:)).^2);
  
  end

  
  
  
end

