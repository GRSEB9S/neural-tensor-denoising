function [X,Zm,Am] = tensorADMM(Y, varargin)
   %[ X,Zm,Am ] = tensorADMM( Y,xnat ) finds underlying rate X for spiking
   %data tensor Y. See Tomioka, Hayashi, Kashima (2010)
   %  Y        4th order tensor, neuron x time x condition x trial. Use NaNs
   %              if condition 1 saw 30 trials but condition 2 saw 28.
   % xnat      n x t x c ground truth for X if available
   %
   % X         n x t x c underlying rate for Y.
   % Zm        3 x 1 cell array for aux variable Z
   % Am        3 x 1 cell array for lagrange multiplier A
   % add param for log(1+exp(x)) to account for smoothness
   
   %% parse inputs
   params = inputParser;
   params.addParamValue('ef', 1, @isscalar); % eta factor;
   params.addParamValue('bf', 1, @isscalar); % beta
   params.addParamValue('gf', 1, @isscalar); % gamma
   params.addParamValue('s',  1, @isscalar); % smoothing param in f(x)
   params.addParamValue('noiseTerm','poiss', @(x) ismember(x,{'poiss','gauss'}));
   params.addParamValue('sig', 1, @isscalar); % sigma for guassian
   params.addParamValue('maxIter', 20, @isscalar);
   params.addParamValue('xnat', []);
   params.parse(varargin{:});
   
   %% copy from params
   ef = params.Results.ef;
   bf = params.Results.bf;
   gf = params.Results.gf;
   s  = params.Results.s;
   noiseTerm = params.Results.noiseTerm;
   sig = params.Results.sig;
   maxIter = params.Results.maxIter;
   xnat = params.Results.xnat;
   
%% ADMM
   method = 1; % which newton method to implement
   eta = 1/(ef*nanstd(Y(:))); % admm parameter, choose factor in front of std(y), see Tomioka (2010) // how to choose? mean(y) instead of std(y) since poisson?
   beta = 1/(bf*nanstd(Y(:))); % total variation parameter // how to choose?
   gamma = gf*ones(3,1); % sum of nuclear norms parameter // how to choose?
   
   [n,t,c,tr] = size(Y);
   
   sz = size(nanmean(Y,4)); % [n t c]
   N = prod(sz); % n*t*c
   k = length(sz); % 3
   ind = (1:N)'; % unused
   [I,J,K] = ind2sub(sz, ind); % unused
   I = {I,J,K}; % unused
   
   prm{1} = [1 2 3]; % permutations
   prm{2} = [2 3 1];
   prm{3} = [3 1 2];
   
   iprm{1} = [1 2 3]; %inverse permutations
   iprm{2} = [3 1 2];
   iprm{3} = [2 3 1];

   %% initialize X,Z,A
   X = nanmean(Y,4);
   X = log(exp(X+0.5) - 1); % 0.5 to keep from exploding to -inf? // possibly better initialization than mean(Y)

   Z = 0 + 0.02*randn(sz); % ADMM aux variable // how to initialize?
   A = 0 + 0.02*randn(sz); % lagrange mult // how to initialize?
   
   Xm = cell(1,k); % unused
   Zm = cell(1,k);
   Am = cell(1,k);
     
   z = cell(k,1); % unused
   a = cell(k,1); % unused
   
   % matricized Z and A (what are actually used in ADMM)
   for kk = 1:k
      Zm{kk} = permute(Z,prm{kk});
      Zm{kk} = Zm{kk}(:,:);
      z{kk} = Zm{kk}(:);
      
      Am{kk} = permute(A,prm{kk});
      Am{kk} = Am{kk}(:,:);
      a{kk} = Am{kk}(:);
   end
     
   % ignore
   %B = zeros(sz); % unused
   %ind = sub2ind(sz, I{:}); % unused
   %B(ind) = y; % unused
   
   %% fr function, log(1+exp(x)) -- consider different functions later
   %s = 1; % can we optimize over s?
   f   = @(x,s) s*log(1+exp(x/s));
   df  = @(x,s) 1./(1+exp(-x/s)); % derivative 
   d2f = @(x,s) exp(x/s)./(s*(1+exp(x/s)).^2); % 2nd derivative

   %%
   %f   = @(x,s) x;
   %df  = @(x,s) ones(size(x));
   %d2f = @(x,s) zeros(size(x));
   
   %% admm

   iter = 0;
   xchange = Inf;
   while xchange > 1e-4 && iter < maxIter % need to implement proper stopping criteria
      stopping = Inf; % stopping criteria for newton
      X_old = X;
      %% x step, newtons method
      switch method
         case 1
            %% method 1 - solve with one big newtons method
            while stopping/norm(X(:),'fro') > 1e-6
               g = zeros(sz); % grad
               h = zeros(sz); % hess

               %% likelihood term
               switch noiseTerm
                  case 'poiss'
                     for rr = 1:tr
                        % old method
                        %yTr = Y(:,:,:,rr);
                        %if sum(isnan(yTr(:))) == 0 
                        %   g = g - Y(:,:,:,rr).*df(X,s)./f(X,s) + df(X,s);
                        %   h = h - Y(:,:,:,rr).*(d2f(X,s).*f(X,s) - df(X,s).^2)./f(X,s).^2 + d2f(X,s);
                        %end
                        % new non shitty method
                        yTr = Y(:,:,:,rr);
                        
                        g_ = -yTr.*df(X,s)./f(X,s) + df(X,s);
                        g_(isnan(g_)) = 0; 
                        g = g + g_;
                        
                        h_ = -yTr.*(d2f(X,s).*f(X,s) - df(X,s).^2)./f(X,s).^2 + d2f(X,s);
                        h_(isnan(h_)) = 0;
                        h = h + h_;
                     end
                  case 'gauss'
                     for rr = 1:tr
                        yTr = Y(:,:,:,rr);
                        if sum(isnan(yTr(:))) == 0
                           % oops, g + 1/sig??....
                           g = g - 1/sig^2*(Y(:,:,:,rr) - f(X,s)).*df(X,s);
                           h = h - 1/sig^2*(Y(:,:,:,rr).*d2f(X,s) - df(X,s).^2 - f(X,s).*d2f(X,s));
                        end
                     end
               end

               %% admm terms
               for kk = 1:k
                  g = g + eta*fold(Am{kk},kk) + eta*(X - fold(Zm{kk},kk)); % add eta to first term or no? this just scales lagrange multiplier
                  h = h + eta*ones(sz);
               end
               
               %% total var term
               dia = 4*ones(t,1); dia(1) = 2; dia(end) = 2; % make a block
               off = -2*ones(t-1,1);
               hblock = diag(dia) + diag(off,1) + diag(off,-1);
               g = g + beta*fold(hblock*flatten(X,2),2); % clever way to get grad, i think
               
               % put t first
               g = permute(g,prm{2});
               h = permute(h,prm{2});
               g = g(:);
               h = h(:); % don't diag(h), just keep track of diag and offdiags
               
               hl = repmat([off; 0],n*c,1); % lower diag
               hu = repmat([off; 0],n*c,1); % upper diag
               hl = beta*hl(1:end-1);       % trim last 0
               hu = beta*hu(1:end-1);       % trim last 0
               hd = h + beta*repmat(dia,n*c,1); % diag of hessian
               
               tn = 1; % step size // implement line search later? after a few admm iterations the x step needs only 1 iteration so probably no need
               xnt = -TDMAsolver([0; hl],hd,[hu; 0],g); % -inv(h)*g // O(n)
               xnt = reshape(xnt,[t c n]); % put back into nxtxc
               xnt = permute(xnt,iprm{2});
               X = X + tn*xnt; % newton step
               stopping = abs(g(:)'*xnt(:)); %stopping criteria
               fprintf('newton <g,xnt>:\t%i\n', stopping);
            end
         case 2
            %% method 2 - do newtons method for each individual tx1 vector (not yet implemented), possibly faster
 
      end
         
      %% z step
      for kk = 1:k
         [U,S,V] = svd(flatten(X,kk) + Am{kk},'econ');
         Zm{kk} = U*max(S-gamma(kk)/eta,0)*V';
      end
      
      %% lagrange mult step
      for kk = 1:k
         Am{kk} = Am{kk} + (flatten(X,kk) - Zm{kk});
      end
      
      xchange = norm(X(:) - X_old(:));
      if ~isempty(xnat)
         xdiff =   norm(X(:) - xnat(:)); else
         xdiff = [];
      end
      iter = iter+1;
      fprintf('iter:\t %i \t x_change:\t %i \n',iter, xchange);
      
   end

   
%% subs  
   function Xout = flatten(Xin,mode)
      % flattens tensor into matrix
      Xin = permute(Xin,prm{mode});
      Xout = reshape(Xin,[sz(mode) prod(sz)/sz(mode)]);
   end

   function Xout = fold(Xin,mode)
      % folds matrix into tensor
      Xin = reshape(Xin,sz(prm{mode}));
      Xout = permute(Xin,iprm{mode});
   end   
   
   function x = TDMAsolver(a,b,c,d)
      % for inverting tri-diagonal hessian
      %a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
      l = length(d);

      % Modify the first-row coefficients
      c(1) = c(1) / b(1);    % Division by zero risk.
      d(1) = d(1) / b(1);   

      for i = 2:l-1
          temp = b(i) - a(i) * c(i-1);
          c(i) = c(i) / temp;
          d(i) = (d(i) - a(i) * d(i-1))/temp;
      end

      d(l) = (d(l) - a(l) * d(l-1))/( b(l) - a(l) * c(l-1));

      % Now back substitute.
      x(l) = d(l);
      for i = l-1:-1:1
          x(i) = d(i) - c(i) * x(i + 1);
      end
   end

end



   
%% trash
%%    
%    IT = 1:N;
%    IT = reshape(IT,[n t c]);
%    
%    it{1} = permute(IT,[1 2 3]);
%    it{2} = permute(IT,[2 3 1]);
%    it{3} = permute(IT,[3 1 2]);
%    
%    it{1} = it{1}(:);
%    it{2} = it{2}(:);
%    it{3} = it{3}(:);
   

%% P
%    I = eye(N);
%    P = cell(k,1);
%    for kk = 1:k
%       P{kk} = I(it{kk},:);
%    end

%% newton method without TV
% while stopping/norm(X(:),'fro') > 1e-6
%    % grad, hess
%    g = zeros(sz);
%    h = zeros(sz);
% 
%    % likelihood
%    for rr = 1:tr
%       yTr = Y(:,:,:,rr);
%       if sum(isnan(yTr(:))) == 0 
%          g = g - Y(:,:,:,rr).*df(X)./f(X) + df(X);
%          h = h - Y(:,:,:,rr).*(d2f(X).*f(X) - df(X).^2)./f(X).^2 + d2f(X);
%       end
%    end
% 
%    % mm terms
%    for kk = 1:k
%       g = g + eta*fold(Am{kk},kk) + eta*(X - fold(Zm{kk},kk)); % add eta to first term or no? this just scales lagrange multiplier
%       h = h + eta*ones(sz);
%    end
% 
%    tn = 1;
%    xnt = -g./h; %newton step
%    X = X + tn*xnt;
%    stopping = abs(g(:)'*xnt(:));
%    %disp(stopping);
% end   
   




