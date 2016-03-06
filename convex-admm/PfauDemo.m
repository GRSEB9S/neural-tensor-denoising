%% Set model parameters
T = 5000; % number of time steps
k = 8; % number of latent states
n = 200; % number of neurons
q = 0; % number of history time steps used (here we have no local connectivity)

A = randn(k);
[v,d] = eig(A);
A = real(v*diag(.99*(cos(angle(diag(d))/20) + 1i*sin(angle(diag(d))/20)))/v); % make the dynamics very slow
[u,s,v] = svd(A);
A = u*min(.99,s)*v'; % Stabilize the dynamics

Q = eye(k);
R = chol(Q)';
C = randn(n,k)/5;
d = -3+randn(n,1); % baseline firing rate
D = 3*randn(n,q*n).*(rand(n,q*n) < 0.05); % 95% sparse connectivity matrix

f = @(x) log(1+exp(x));

%% Create model data
x = zeros(k,T+q); % latent state
e = R*rand(k,T+q)/2;%R*randn(k,T+q);
x(:,1) = e(:,1);
L = zeros(n,T+q);
X = zeros(n,T+q);
for t = 2:T+q
    x(:,t) = A*x(:,t-1) + e(:,t);
    L(:,t) = C*x(:,t) + d;
    if t > q
        L(:,t) = L(:,t) + D*vec(X(:,t-q:t-1));
    end
    X(:,t) = poissrnd(f(L(:,t)));
end

figure
imagesc(X), colorbar, title('Spike Train')

%% Run nuclear norm minimization
opts.q = q;
opts.verbose = 1;
opts.lambda = .001;
opts.gamma = 0;
opts.nlin = 'logexp'; % if f = @exp, replace with 'exp'

disp('Run Nuc Norm Min')

%Y = nucnrmmin(X,opts);
Xout = nucnrmminSimple(X,opts);
figure
plot(svd(Y(sum(X,2)~=0,:)))
title('Spectrum')