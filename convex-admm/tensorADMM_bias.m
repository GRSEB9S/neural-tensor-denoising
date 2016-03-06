function [X, B] = tensorADMM_bias(Y, varargin)
%X = tensorADMM( Y, varargin ) finds underlying rate X for spiking
%data tensor Y with a convex optimization approach. The problem is solved
%using the alternating direction method of multipliers (ADMM). The loss
%function is
%
%  L = -sum log p(y|x) + beta*smoothness_term + gamma*tensor_rank_term
%
% -sum log p(y|x)  -- log likelihood term (poisson) % each neuron has a bias
% smoothness_term  -- term to ensure smoothness of signals (useful for low
%                     FR neurons
% tensor_rank_term -- term to penalize multilinear rank
%
% See Tomioka, Hayashi, Kashima (2010)
%
% INPUTS 
% Y            4th order tensor, neuron x time x condition x trial. Use NaNs
%                 if condition 1 saw 30 trials but condition 2 saw 28, etc
% 'ef'         eta parameter in augmented lagrangian.
% 'bf'         beta regularization parameter for total variation term.
% 'gf'         gamma regularization parameter for trace norm terms, lower
%                 for smoother signals
% 's'          how smooth to make the log(1 + exp()) function. default 1,
%              smaller values for a sharper threshold. Too sharp and
%              Newton's method won't work (careful).
% 'maxIter'    max number of iterations. 
%
% OUTPUTS
% X         n x t x c underlying rate for Y.
% B         n bias terms for firing rates
%
% NOTE: this code is fully functioning but is still a work in progress.
% please contact me if you intend to use this and i can provide an updated
% version -- jsseely@gmail.com
%
% Jeffrey Seely, Columbia University
% Bias term added by Josh Merel

%% parse inputs
   params = inputParser;
   params.addParamValue('ef', 1, @isscalar); % eta factor;
   params.addParamValue('bf', 1, @isscalar); % beta
   params.addParamValue('gf', 1, @isscalar); % gamma
   params.addParamValue('s',  1, @isscalar); % smoothing param in f(x)
   params.addParamValue('noiseTerm','poiss', @(x) ismember(x,{'poiss','gauss'})); % for later
   params.addParamValue('sig', 1, @isscalar); % sigma for guassian (for later)
   params.addParamValue('maxIter', 20, @isscalar);
   params.parse(varargin{:});
   
%% copy from params
   ef = params.Results.ef;
   bf = params.Results.bf;
   gf = params.Results.gf;
   s  = params.Results.s;
   noiseTerm = params.Results.noiseTerm;
   sig = params.Results.sig;
   maxIter = params.Results.maxIter;
   
%% regularization params, etc
   method = 1; % which newton method to implement -- only method 1 implemented so far
   eta = 1/(ef*nanstd(Y(:))); % admm parameter, choose factor in front of std(y), see Tomioka (2010)
   beta = 1/(bf*nanstd(Y(:))); % total variation parameter
   gamma = gf*ones(3,1); % sum of nuclear norms parameter
   
   [n,t,c,tr] = size(Y);
   
   sz = size(nanmean(Y,4)); % [n t c]
   N = prod(sz); % n*t*c
   k = length(sz); % 3
   
   prm{1} = [1 2 3]; % permutations
   prm{2} = [2 3 1];
   prm{3} = [3 1 2];
   
   iprm{1} = [1 2 3]; %inverse permutations
   iprm{2} = [3 1 2];
   iprm{3} = [2 3 1];

%% initialize X,Z,A
   X = nanmean(Y,4);
   X = log(exp(X+0.5) - 1); % 0.5 to keep from exploding to -inf
   B = zeros(n,1); 
   
   Z = 0 + 0.02*randn(sz); % ADMM aux variable
   A = 0 + 0.02*randn(sz); % lagrange mult
   
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
     
%% fr function, log(1+exp(x))
%    f   = @(x,s) s*log(1+exp(x/s));
%    df  = @(x,s) 1./(1+exp(-x/s)); % derivative 
%    d2f = @(x,s) exp(x/s)./(s*(1+exp(x/s)).^2); % 2nd derivative

   f   = @(x,bias,s) s*log(1+exp(bsxfun(@plus, bias, x)/s));
   df  = @(x,bias,s) 1./(1+exp(-bsxfun(@plus, bias, x)/s)); % derivative 
   d2f = @(x,bias,s) exp(bsxfun(@plus, bias, x)/s)./(s*(1+exp(bsxfun(@plus, bias, x)/s)).^2); % 2nd derivative
   
%% 
   %f   = @(x,s) x;
   %df  = @(x,s) ones(size(x));
   %d2f = @(x,s) zeros(size(x));
   
%% admm
   iter = 0;
   xchange = Inf;
   while xchange > 1e-4 && iter < maxIter 
      stopping = Inf; % stopping criteria for newton
      X_old = X;
        
      %% x step, newtons method
      switch method
         case 1
            %% method 1 - solve with one big newtons method
            while stopping/norm(X(:),'fro') > 1e-2
               g = zeros(sz); % grad
               h = zeros(sz); % hess
                
               %% likelihood term
               switch noiseTerm
                  case 'poiss'
                     for rr = 1:tr
                        yTr = Y(:,:,:,rr);
                        if sum(isnan(yTr(:))) == 0 
                           g = g - Y(:,:,:,rr).*df(X,B,s)./f(X,B,s) + df(X,B,s);
                           h = h - Y(:,:,:,rr).*(d2f(X,B,s).*f(X,B,s) - df(X,B,s).^2)./f(X,B,s).^2 + d2f(X,B,s);
                        end
                     end
                  case 'gauss'
                     % for later... sorry guass
               end
               %bias gradient and hessian
               g_b = sum(sum(g,2),3);
               h_b = h; %sum(sum(h,2),3);

               %% admm terms
               for kk = 1:k
                  g = g + eta*fold(Am{kk},kk) + eta*(X - fold(Zm{kk},kk)); % add eta to first term or no? this just scales lagrange multiplier
                  h = h + eta*ones(sz);
               end
               
               %% total var term
               dia = 4*ones(t,1); dia(1) = 2; dia(end) = 2; % make a block
               off = -2*ones(t-1,1);
               hblock = diag(dia) + diag(off,1) + diag(off,-1);
               g = g + beta*fold(hblock*flatten(X,2),2); % nice way to get grad
               
               % put t first *****
               g = permute(g,prm{2});
               h = permute(h,prm{2});
               g = g(:);
               h = h(:); % don't diag(h), just keep track of diag and offdiags
               
               hl = repmat([off; 0],n*c,1); % lower diag
               hu = repmat([off; 0],n*c,1); % upper diag
               hl = beta*hl(1:end-1);       % trim last 0
               hu = beta*hu(1:end-1);       % trim last 0
               hd = h + beta*repmat(dia,n*c,1); % diag of hessian
               
               tn = 1; % step size // implement line search later // after a few admm iterations the x step needs only 1 iteration so probably no need

               %%%%% TDMA solver (for just H_xx) %%%%%%%
%                xnt = -TDMAsolver([0; hl],hd,[hu; 0],g); % -inv(h)*g // O(n)
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               %added by Josh
               %%%%%%% sparse matrix solver for H_bx,H_xx,H_bb jointly %%%%%%%
               if ~exist('tcnInd','var')
%                    ntcInd = reshape(1:(n*t*c),n,t,c);
                    tcnInd = reshape(1:(n*t*c),t,c,n);
               end
               if ~exist('b_cross','var')
                   b_cross = repmat((1:n)',1,t,c);
                   b_cross = permute(b_cross,prm{2});
               end
               %try big sparse Hessian for B and X jointly
               H_xx = spdiags([[hl;0] hd [0;hu]], -1:1, length(hd), length(hd));
               %h_b  is in  n,T,c form and want to fill in n x n*T*c matrix
               H_bb = spdiags(sum(sum(h_b,2),3),0,n,n);
               h_b = permute(h_b,prm{2});
               H_bx = sparse(b_cross(:),tcnInd(:),h_b(:));
               xnt = - [H_xx H_bx'; H_bx H_bb]\[g; g_b];
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               B_upd = xnt((end-n+1):end); 
               xnt = reshape(xnt(1:(end-n)),[t c n]); % put back into nxtxc
               xnt = permute(xnt,iprm{2});
               X = X + tn*xnt; % newton step
               B = B + tn*B_upd;
               
               % for quite large n,t,c, perhaps we need to separately do
               % Newton update for each neuron?
               
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
      iter = iter+1;
      fprintf('iter:\t %i \t x_change:\t %i \n',iter, xchange);
      
   end

   
%% subfunctions
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

