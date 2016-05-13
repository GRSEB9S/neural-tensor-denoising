%% Quantifying condition-trial tradeoff.
% simulate one neuron for C conditions
clear; clc; close all;
%% 
C = 50; % conditions
R = 2;   % trials
T = 100; % time points

%% Truth. Specify condition space.
t_axis = linspace(0,2*pi,T);
freqs = linspace(1,2,C); % frequencies -- range from _ cycles to _ cycles
amps = randn(1,C);

m = repmat(amps', [1 T]).*sin(bsxfun(@plus, freqs'*t_axis, 2*pi*rand(C,1))); % C x T signals (true means)
m_1 = m(1,:); % focus on condition 1
S_c = cov(m'); 

close all; figure;
plot(m') % plot all condition means
figure;
plot(svds(S_c,C)) % rank across conditions

disp(rank(S_c,0.01)) % disp approx rank

%% TxT Cov: (Same for all conditions?)
% Give S_t no structure
%S_t = randn(T);
%S_t = S_t'*S_t/(T);

% Give S_t structure
S_t = 1*cov(sin(bsxfun(@plus, linspace(0.1,10,C)'*t_axis, 2*pi*rand(C,1)))); % give S_t something reasonable

close all;
imagesc(S_t);
figure;
plot(svds(S_t,T))

%% Prior:
m_0 = zeros(T,1);
S_0 = eye(T);

%% data
X = [];
for cc = 1:C
    X(:,:,cc) = mvnrnd(m(cc,:), S_t, R*C);
end

%% Scenario 1: one condition, R trials
X_1 = X(1:R,:,1);
close all; figure; hold all;
plot(X_1','k.')
plot(m_1,'r','linewidth',3)
plot(mean(X_1),'k','linewidth',3)

err_1 = norm(mean(X_1) - m_1)^2;

%% Scenario 2: one condition, RC trials
X_2 = X(:,:,1);
close all; figure; hold all; 
plot(X_2','k.')
plot(m_1,'r','linewidth',3)
plot(mean(X_2),'k','linewidth',3)

err_2 = norm(mean(X_2) - m_1)^2;

%% Scenario 3: C conditions, R trials
    X_3 = X(1:R,:,:);

    %% SVD.
    X_3_mean = squeeze(mean(X_3,1));
    [u,s,v] = svd(X_3_mean);
    X_rec = [];
    for kk = 1:min(T,C)
        X_rec{kk} = u(:,1:kk)*s(1:kk,1:kk)*v(:,1:kk)';
    end

    %% plot 
    close all; figure; hold all;
    plot(X_3(:,:,1)','k.');
    plot(m_1,'r','linewidth',3);
    plot(X_3_mean(:,1),'k','linewidth',2);

    disp(rank(S_c,0.01)) % disp approx rank 
    
    err_3 = [];
    for kk = 1:min(T,C)
        err_3(kk) = norm(X_rec{kk}(:,1) - m_1')^2;
    end
    % plot rank kk reconstruction
    kk = find(err_3 == min(err_3),1);
    plot(X_rec{kk}(:,1),'b','linewidth',3);

    plot(mean(X_2),'k','linewidth',2);
    
    h = figure; hold all;
    plot(err_3,'b')
    plot([1 C], [err_1 err_1],'k'); % upper bound
    plot([1 C], [err_2 err_2],'k'); % lower bound
    
    
%% svd version 2
% X_3_fl = reshape(permute(X_3,[2 1 3]), R*T, C);
% [u,s,v] = svd(X_3_fl, 'econ');
% X_rec = [];
% for kk = 1:C
%     X_rec{kk} = u(:,1:kk)*s(1:kk,1:kk)*v(:,1:kk)';
%     X_rec{kk} = reshape(X_rec{kk}, T,R,C);
%     X_rec{kk} = permute(X_rec{kk}, [2 1 3]);
% end
% 
% %%
% close all; figure; hold all;
% plot(X_3(:,:,1)','k.');
% plot(m_true,'r','linewidth',3);
% plot(X_3_mean(:,1),'k');
% 
% for kk = 1
%     plot(squeeze(mean(X_rec{kk}(:,:,1),1)),'b','linewidth',1);
% end
























