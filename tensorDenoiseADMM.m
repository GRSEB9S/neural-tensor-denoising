function [X, B] = tensorDenoiseADMM(Y, varargin)
% 2 denotes that this is a major update (2015) to the 2014 version.
% goal is to include both 'overlapped' and 'latent' Schatten 1-norm options
% as well as bias term (from Josh).
%X = tensorADMM( Y, varargin ) finds underlying rate X for spiking
%data tensor Y with a convex optimization approach. The problem is solved
%using the alternating direction method of multipliers (ADMM). The loss
%function is
%
%  L = -sum log p(y|f(x)) + beta*smoothness_term + gamma*tensor_rank_term
%
% -sum log p(y|x)  -- log likelihood term (poisson); each neuron has a bias
% smoothness_term  -- term to ensure smoothness of signals (useful for low
%                     FR neurons)
% tensor_rank_term -- term to penalize multilinear rank
%
% See Tomioka, Hayashi, Kashima (2010)
%
% INPUTS 
% Y            4th order tensor, neuron x time x condition x trial. Use NaNs
%                 if condition 1 saw 30 trials but condition 2 saw 28, etc
% 'ef'         eta parameter in augmented lagrangian.
% 'bf'         beta regularization parameter for smoothness term.
% 'lf'         lambda regularization parameter for trace norm terms, lower
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
params.addParameter('ef', 1, @isscalar); % eta
params.addParameter('bf', 1, @isscalar); % beta
params.addParameter('lf', 1, @isscalar); % lambda
params.addParameter('s',  1, @isscalar); % smoothing param in f(x)
params.addParameter('noiseTerm','poiss', @(x) ismember(x,{'poiss','gauss'})); % for later
params.addParameter('smoothingMethod','justX', @(x) ismember(x,{'justX','postNL'}));
params.addParameter('sig', 1, @isscalar); % sigma for guassian (for later)
params.addParameter('maxIter', 20, @isscalar);
params.addParameter('trueX', [], @isnumeric);
params.addParameter('trueB', [], @isnumeric);
params.parse(varargin{:});

%% copy from params
ef = params.Results.ef;
bf = params.Results.bf;
lf = params.Results.lf;
s  = params.Results.s;
noiseTerm = params.Results.noiseTerm;
smoothingMethod = params.Results.smoothingMethod;
sig = params.Results.sig;
maxIter = params.Results.maxIter;
trueX = params.Results.trueX;
trueB = params.Results.trueB;

%% regularization params, etc
method = 1; % which newton method to implement -- only method 1 implemented so far
eta = 1/(ef*std(Y(:),'omitnan')); % admm parameter, choose factor in front of std(y), see Tomioka (2010)
%beta = 1/(bf*std(Y(:),'omitnan')); % smoothing parameter
beta = bf; % smoothing parameter? 
lambda = lf; % sum of nuclear norms parameter

[n,t,c,tr] = size(Y);

sz = size(mean(Y,4,'omitnan')); % [n t c]
N = prod(sz); % n*t*c
k = length(sz); % 3

prm{1} = [1 2 3]; % permutations
prm{2} = [2 3 1];
prm{3} = [3 1 2];

% see ipermute function?
iprm{1} = [1 2 3]; %inverse permutations
iprm{2} = [3 1 2];
iprm{3} = [2 3 1];

%% initialize X,Z,A
X = nanmean(Y,4); % rates
X = log(exp(X+0.5) - 1); % 0.5 to keep from exploding to -inf
B = zeros(n,1); % bias terms

Z = 0 + 0.02*randn(sz); % ADMM aux variable
A = 0 + 0.02*randn(sz); % lagrange mult

Zm = cell(1,k);
Am = cell(1,k);

% matricized Z and A (what are actually used in ADMM)
for kk = 1:k
   Zm{kk} = permute(Z,prm{kk});
   Zm{kk} = Zm{kk}(:,:);

   Am{kk} = permute(A,prm{kk});
   Am{kk} = Am{kk}(:,:);
end
     
%% fr function, log(1+exp(x + b))
% f   = @(x,s) s*log(1+exp(x/s));
% df  = @(x,s) 1./(1+exp(-x/s)); % derivative 
% d2f = @(x,s) exp(x/s)./(s*(1+exp(x/s)).^2); % 2nd derivative

f   = @(x,b,s) s*log(1+exp(bsxfun(@plus, b, x)/s));
df  = @(x,b,s) 1./(1+exp(-bsxfun(@plus, b, x)/s)); % derivative w.r.t. x or w.r.t. b
d2f = @(x,b,s) exp(bsxfun(@plus, b, x)/s)./(s*(1+exp(bsxfun(@plus, b, x)/s)).^2); % 2nd derivative w.r.t. x or w.r.t b

%f   = @(x,s) x;
%df  = @(x,s) ones(size(x));
%d2f = @(x,s) zeros(size(x));

%% plot init
rows = 5;
cols = 5;
h1 = zeros(rows*cols,1);
h2 = zeros(rows*cols,1);
h3 = zeros(rows*cols-1,1);
plotInit
plotInline(0)

%% admm
iter = 0;
xchange = Inf;
while xchange > 1e-1 && iter < maxIter % better stopping criteria? 1e-4 is too low.
   stopping = Inf; % stopping criteria for newton
   X_old = X;
   B_old = B;

   %% x step, newtons method
   switch method
   case 1
   %% method 1 - sum of nuclear norms / 'overlapped' schatten 1 norm
      while stopping/norm([X(:); B],'fro') > 1e-1
         g_x = zeros(sz); % gradient terms w.r.t. x
         h_xx = zeros(sz); % hessian diagonal terms w.r.t. xx

         %% likelihood term, g_x and h_xx
         switch noiseTerm
         case 'poiss'
            for rr = 1:tr
               y_tr = Y(:,:,:,rr);
               g_tr = -y_tr.*df(X,B,s)./f(X,B,s) + df(X,B,s);
               g_tr(isnan(g_tr)) = 0;
               g_x = g_x + g_tr;
               
               h_tr = -y_tr.*(d2f(X,B,s).*f(X,B,s) - df(X,B,s).^2)./f(X,B,s).^2 + d2f(X,B,s);
               h_tr(isnan(h_tr)) = 0;
               h_xx = h_xx + h_tr;
            end
         case 'gauss'
            for rr = 1:tr
               y_tr = Y(:,:,:,rr);
               % 1/sig or 1/sig^2?
               % double check ---
               g_tr = -1/sig^2*(y_tr - f(X,B,s)).*df(X,B,s);
               g_tr(isnan(g_tr)) = 0;
               g_x = g_x + g_tr;
               
               h_tr = -1/sig^2*(y_tr.*d2f(X,B,s) - df(X,B,s).^2 - f(X,B,s).*d2f(X,B,s));
               h_tr(isnan(h_tr)) = 0;
               h_xx = h_xx + h_tr;
            end
         end

         %% smooth term, g_x and h_xx
         switch smoothingMethod
         case 'postNL'
            % double check -- is convex??
            % c_mat contains integer coefficients useful in defining grad
            % and hess for smoothness term
            c_diag = 4*ones(t,1); c_diag(1) = 2; c_diag(end) = 2; % diagonal term for c_mat
            c_off = -2*ones(t-1,1); % off diagonal terms for c_mat
            c_mat = diag(c_diag) + diag(c_off,1) + diag(c_off,-1); % create block to get gradient

            g_xs = zeros(sz);
            h_xxs = zeros(sz); % just the diag component
            h_off = zeros(sz - [0 1 0]); % off-diag component of h_xxs
            for nn = 1:n
            for cc = 1:c
               g_xs(nn,:,cc) = diag(df(X(nn,:,cc),B(nn),s))*c_mat*f(X(nn,:,cc),B(nn),s)';
               h_xxs(nn,:,cc) = diag(d2f(X(nn,:,cc),B(nn),s))*c_mat*f(X(nn,:,cc),B(nn),s)' + c_diag.*df(X(nn,:,cc),B(nn),s)'.^2; % diag component of h_xxs
               h_off(nn,:,cc) = c_off.*df(X(nn,1:end-1,cc),B(nn),s)'.*df(X(nn,2:end,cc),B(nn),s)'; % off diag component of h_xxs
            end
            end

            % useful
            h_upper = beta*cat(2, zeros(n,1,c), h_off);
            h_lower = beta*cat(2, h_off, zeros(n,1,c));

            % update g_x and h_xx
            g_x = g_x + beta*g_xs;
            h_xx = h_xx + beta*h_xxs; % diag component only
         
            % LL term and smooth term for g_b, h_bb, h_bx
            g_b = sum(sum(g_x,2),3);
            h_bb = sum(sum(h_xx,2),3);
            h_bx = h_xx + h_lower + h_upper;
         case 'justX'
            % LL term for g_b, h_bb, h_bx
            g_b = sum(sum(g_x,2),3);
            h_bb = sum(sum(h_xx,2),3);
            h_bx = h_xx;
            
            c_diag = 4*ones(t,1); c_diag(1) = 2; c_diag(end) = 2; % diagonal term for c_mat
            c_off = -2*ones(t-1,1); % off diagonal terms for c_mat
            c_mat = diag(c_diag) + diag(c_off,1) + diag(c_off,-1);
            g_xs = fold(c_mat*flatten(X,2),2);

            h_xxs = repmat(c_diag', n, 1, c); 
            h_upper = -beta*2*ones(sz);
            h_lower = -beta*2*ones(sz);
            h_upper(:,1,:) = 0;
            h_lower(:,end,:) = 0;
            
            g_x = g_x + beta*g_xs;
            h_xx = h_xx + beta*h_xxs; % diag entries only
         end
                  
         %% LL term and smooth term; g_b, h_bb, h_bx
         %g_b = sum(sum(g_x,2),3);
         %h_bb = sum(sum(h_xx,2),3);
         %h_bx = h_xx + h_lower + h_upper;
        
         %% admm terms; g_x, h_xx
         for kk = 1:k
            % no need to fold/refold??? just 3*eta*A??
            g_x = g_x + eta*fold(Am{kk},kk) + eta*(X - fold(Zm{kk},kk)); % add eta to first term or no? this just scales lagrange multiplier // 
            h_xx = h_xx + eta*ones(sz);
         end
         
         %% solve; newton step
         x_upd = zeros(sz);
         b_upd = zeros(n,1);
         % solve neuron by neuron
         for nn = 1:n
            H = spdiags([h_lower(nn,:)' h_xx(nn,:)' h_upper(nn,:)'], -1:1, t*c, t*c); % correct?
            Hbx = sparse(h_bx(nn,:));
            Hbb = sparse(h_bb(nn));
            upd = -[H Hbx'; Hbx Hbb]\[g_x(nn,:)'; g_b(nn)];
            x_upd(nn,:,:) = reshape(upd(1:end-1), [1 t c]);
            b_upd(nn) = upd(end);
         end
         
         % newton step
         tn = 1;
         X = X + tn*x_upd;
         B = B + tn*b_upd;
         
         plotInline(1)
         stopping = abs([g_x(:); g_b(:)]'*[x_upd(:); b_upd]); %stopping criteria
         fprintf('newton <g,xnt>:\t%i\n', stopping);
         
      end
   case 2
   %% method 2 - latent schatten-1 norm

   end

   %% z step
   for kk = 1:k
      [U,S,V] = svd(flatten(X,kk) + Am{kk},'econ'); % - or + ???? differences in tomioka papers...
      Zm{kk} = U*max(S-lambda/eta,0)*V';
   end

   %% lagrange mult step
   for kk = 1:k
      Am{kk} = Am{kk} + (flatten(X,kk) - Zm{kk}); % use eta here?
   end

   %% wrapup
   plotInline(1)
   xchange = norm([X(:);B] - [X_old(:);B_old]);
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

function plotInit
   % initialize plotting
   close all
   subs = rows*cols-1;
   figure(1);
   for sub_ind = 1:subs
      [nind, cind] = ind2sub([rows cols],sub_ind);
      subplot(rows,cols,sub_ind);
      xlim([1 t]);
      hold on;
      h1(sub_ind) = plot(mean(Y(nind,:,cind,:),4,'omitnan'),'k','linewidth',2);
      if ~isempty(trueX)
         plot(trueX(nind,:,cind),'b','linewidth',2);
      end
   end
   subplot(rows,cols,rows*cols);
   xlim([1 n]);
   hold on;
   if ~isempty(trueB)
      h1(end) = plot(trueB,'k','linewidth',2); 
   end
end

function plotInline(delh)
   % plot x and b estimates
   subs = rows*cols-1;
   figure(1);
   for spind = 1:subs
      [nind, cind] = ind2sub([rows cols],spind);
      subplot(rows,cols,spind);
      if delh
         delete(h2(spind));
         delete(h3(spind));
      end
      h2(spind) = plot(f(X(nind,:,cind),B(nind),s),'r','linewidth',1.5);
      h3(spind) = plot(X(nind,:,cind),'g','linewidth',1.5);
      drawnow
   end
   subplot(rows,cols,rows*cols);
   if delh
      delete(h2(end));
   end
   h2(end) = plot(B,'r');
end
   
end

% wolfram alpha example for checking derivatives:
% x:=x_t-1, y:=x_t, z:=x_t+1
%
% d/dy d/db [(log(1+exp(z+b)) - log(1+exp(y+b)))^2 + (log(1+exp(y+b)) - log(1+exp(x+b)))^2]
%





