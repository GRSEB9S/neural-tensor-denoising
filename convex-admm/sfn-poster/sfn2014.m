%% SFN 2014 -- analysis script
%% load RC
clear
RCpath = 'Datasets/Motor/RC,2009-09-18,1-2,good-ss';
Spath = 'Datasets/Motor/Lara-20141105/Sstructs/';

%% pre-process S structs (Lara)
bw = 10;
sd = 10;

[T1,A1,S1] = preS(Spath, bw, sd);
%[~,~] = preRC(RCpath, bw, sd);

bw = 25;
sd = 25;

[T2,A2,S2] = preS(Spath, bw, sd); % -- REDO SO BW WORKS (now only takes factors of 1300);
%[~,~] = preRC(RCpath, bw, sd);

%% cleanup
clearvars -except T1 A1 S1 T2 A2 S2 RCpath Spath

      %% addm sandbox
%       s = .05;
%       f   = @(x,s) s*log(1+exp(x/s));
% 
%       %% test tensor_admm
%       T_admm = tensorADMM(T2, 'ef', 0.01, 'bf', .05, 'gf', 0, 's', s, 'maxIter',5);
% 
%       %% plot results
%       %plotAgrid(tensorOpt(f(T_admm,s), 70, 1));
%       plotAgrid(f(T_admm,s))
%       plotAgrid(nanmean(T2,4));
%       
%       %% 
%       T0_admm = tensorADMM(T_tr{1}, 'ef', 0.01, 'bf', 0.1, 'gf', 0, 's', s, 'maxIter', 20);
%       %%
%       plotAgrid(f(T0_admm,s))
%       plotAgrid(nanmean(T_tr{1},4));
      
      
%% main analysis - Lara data
bw = 10; % using A1, hence bw = 10;

[n,t_t,c,r] = size(T2);
t_a = size(A1,2);

for nn = 1:n
   rperm = randperm(r);
   T2(nn,:,:,:) = T2(nn,:,:,rperm);
   A1(nn,:,:,:) = A1(nn,:,:,rperm);
   T1(nn,:,:,:) = T1(nn,:,:,rperm);
end

% params
P = 10; % number of params to search over
T2param = linspace(0,50,P);  %old: 1 to 100. new: 0 to 50
A1param = round(linspace(30,100,P));
A0param = linspace(10,100,P);
s = 0.5; % old: 0.05 -- way too small... new: 0.5?
f   = @(x,s) s*log(1+exp(x/s));
f_  = @(x) x.*(x>0); % not used
ef = 0.01;
bf = 0.05;
maxIter = 20;

% cross val
K = 5;
b = round((0:K)*r/K);

for kk = 1:K

   % indices
   b_tr = 1:r;
   b_val = (b(kk)+1):(b(kk+1));
   b_tr(b_val) = [];
   
   % data - cells usage is temporary
   T2_tr{kk}  = T2(:,:,:,b_tr);
   T1_tr{kk}  = T1(:,:,:,b_tr);
   A1_tr{kk}  = A1(:,:,:,b_tr);
   T2_val{kk} = T2(:,:,:,b_val);
   T1_val{kk} = T1(:,:,:,b_val);
   A1_val{kk} = A1(:,:,:,b_val);
   
   % null model
   %T0_admm{kk} = tensorADMM(T_tr{kk}, 'ef', ef, 'bf', bf, 'gf', 0, 's', s, 'maxiter', maxIter);
   
   2str(pp)]);
      % models
      T2_admm{pp,kk} = tensorADMM(T2_tr{kk}, 'ef', ef, 'bf', bf, 'gf', T2param(pp), 's', s, 'maxiter', maxIter);
      A1_tens{pp,kk} = double(tensorOpt(nanmean(A1_tr{kk},4), A1param(pp), 1));

      % null models
      [~,A0_tens{pp,kk},~] = preS(Spath, bw, A0param(pp)); % bw should be 10.
      A0_tens{pp,kk} = nanmean(A0_tens{pp,kk},4);
      
      % errors
      T2_admm_rep = repmat(T2_admm{pp,kk},1,1,1,size(T2_val{kk},4));
       
      A1_tens_rep = repmat(A1_tens{pp,kk},1,1,1,size(A1_val{kk},4));
      A1_tens_rep = A1_tens_rep*(bw/1000);
      
      A0_tens_rep= repmat(A0_tens{pp,kk},1,1,1,size(A1_val{kk},4));
      A0_tens_rep = A0_tens_rep*(bw/1000);

      % non-nan indices
      nn_inds_t2 = ~isnan(T2_val{kk});
      nn_inds_t1 = ~isnan(T1_val{kk});
      nn_inds_a1 = ~isnan(A1_val{kk}); % should be the same

      % errors (- log likelihood)
      T2_ll_loop = -log(vec(poisspdf(...
                    T2_val{kk}(nn_inds_t2), f(T2_admm_rep(nn_inds_t2),s) ...
                    )));
      A1_ll_loop = -log(vec(poisspdf(...
                    T1_val{kk}(nn_inds_t1), f(A1_tens_rep(nn_inds_t1),s) ...
                    ))); % using f(,s) on A is not totally correct... oh well fuck it
      A0_ll_loop = -log(vec(poisspdf(...
                    T1_val{kk}(nn_inds_t1), f(A0_tens_rep(nn_inds_t1),s) ...
                    )));
               
      % errors (MSE)
      T2_mse_loop = vec(T2_val{kk}(nn_inds_a1) - f(T2_admm_rep(nn_inds_t2),s)).^2;
      A1_mse_loop = vec(A1_val{kk}(nn_inds_a1) - A1_tens_rep(nn_inds_a1)).^2;
      A0_mse_loop = vec(A1_val{kk}(nn_inds_a1) - A0_tens_rep(nn_inds_a1)).^2;
      
      % means of errors
      T2_ll(pp,kk) = nanmean(T2_ll_loop);
      A1_ll(pp,kk) = nanmean(A1_ll_loop);
      A0_ll(pp,kk) = nanmean(A0_ll_loop);
      T2_mse(pp,kk) = nanmean(T2_mse_loop);
      A1_mse(pp,kk) = nanmean(A1_mse_loop);
      A0_mse(pp,kk) = nanmean(A0_mse_loop);
      
   end
end



%% some plotting
nn = 97;
cc = 6;
figure; hold on;
plot(nanmean(T1(nn,:,cc,:),4),'linewidth',2);
for kk = 1:5
   plot(linspace(1,130,52),52/130*f(T_admm{1,kk}(nn,:,cc),.05));
   plot(f(A_tens{1,kk}(nn,:,cc),10)*(bw/1000))
end

figure; hold on;
plot(nanmean(T2(nn,:,cc,:),4),'linewidth',2);
for kk = 1:5
   plot(f(T_admm{1,kk}(nn,:,cc),0.5));
end

      %% addm sandbox
      s = .5;
      f   = @(x,s) s*log(1+exp(x/s));

      %% test tensor_admm
      T_test = tensorADMM(T2, 'ef', 0.01, 'bf', 0.5, 'gf', 70, 's', s, 'maxIter',50);

      %% plot results
      nn = 97;
      cc = 6;
      s = 0.5;
      figure; hold on;
      plot(nanmean(T2(nn,:,cc,:),4),'linewidth',2);
      plot(f(T_test(nn,:,cc),s));

