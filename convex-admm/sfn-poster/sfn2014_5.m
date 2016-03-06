%% sfn2014_5 -- 5-fold CV on 5 trials only -- test on remaining trials
%% load data
clear
load('sfn_10bw_10sd.mat');
%% prepare arrays
A  = D.L.A;
T  = D.L.T;
T0 = D.L.T0;
bw = D.bw;

%% train/val vs test sets
A_test = A(:,:,:,6:end);
T_test = T(:,:,:,6:end);
T0_test=T0(:,:,:,6:end);

A = A(:,:,:,1:5);
T = T(:,:,:,1:5);
T0=T0(:,:,:,1:5);

[n,t,c,r] = size(T);

%% main analysis
% randomize neuron by neuron
for nn = 1:n
   rperm = randperm(r);
   T(nn,:,:,:) = T(nn,:,:,rperm);
   A(nn,:,:,:) = A(nn,:,:,rperm);
   T0(nn,:,:,:) = T0(nn,:,:,rperm);
end

% params
P = 50;
Aparam = round(linspace(3,150,P));
A0param = linspace(1,10,P);
s = 0.1;
f = @(x,s) s*log(1+exp(x/s));

% cross val
K = 5;
b = round((0:K)*r/K);

%% main loop
for kk = 1:K

   % indices
   b_tr = 1:r;
   b_val = (b(kk)+1):(b(kk+1));
   b_tr(b_val) = [];
   
   % data - cells usage is temporary
   T_tr{kk}  = T(:,:,:,b_tr);
   A_tr{kk}  = A(:,:,:,b_tr);
   T_val{kk} = T(:,:,:,b_val);
   A_val{kk} = A(:,:,:,b_val);
   
   for pp = 1:P
      disp(['Models. k = ' num2str(kk) ', p = ' num2str(pp)])
      % models
      A_tens{pp,kk} = double(tensorOpt(nanmean(A_tr{kk},4), Aparam(pp), 1));

      % null model
      A0_loop = filterGauss(A_tr{kk}, A0param(pp), 2);
      %A0_loop = A0_loop(:,1:bw:end,:,:);
      A0_loop = nanmean(A0_loop,4);
      A0_tens{pp,kk} = A0_loop;
      
   end
end   

%% compute errors
for kk = 1:K
for pp = 1:P
   
   disp(['Errors. k = ' num2str(kk) ', p = ' num2str(pp)])

   bf = bw/1000; % bin factor
   %A_ll_val(pp,kk)   = sfn_ll( T_val{kk}, f(bf*A_tens{pp,kk},s) );
   %A0_ll_val(pp,kk)  = sfn_ll( T_val{kk}, f(bf*A0_tens{pp,kk},s) );
   A_mse_val(pp,kk)  = sfn_mse( A_val{kk}, A_tens{pp,kk} );
   A0_mse_val(pp,kk) = sfn_mse( A_val{kk}, A0_tens{pp,kk} );
   
   %A_ll_tr(pp,kk) = sfn_ll( T_tr{kk}, f(bf*A_tens{pp,kk},s) );
   %A0_ll_tr(pp,kk) = sfn_ll( T_tr{kk}, f(bf*A0_tens{pp,kk},s) );
   A_mse_tr(pp,kk) = sfn_mse( A_tr{kk}, A_tens{pp,kk} );
   A0_mse_tr(pp,kk) = sfn_mse( A_tr{kk}, A0_tens{pp,kk} );
   
end
end

%% model assessment
[a_mse_val,a_i] = min(mean(A_mse_val,2));
[a0_mse_val,a0_i]=min(mean(A0_mse_val,2));

%%
A_tens_best = double(tensorOpt(nanmean(A,4), Aparam(a_i), 1));
A0_tens_best= filterGauss(A, A0param(a0_i), 2);
A0_tens_best= nanmean(A0_tens_best,4);

A_mse_test = sfn_mse( A_test, A_tens_best );
A0_mse_test= sfn_mse( A_test, A0_tens_best);


%%

figure; hold on;
plot(mean(A_mse_val,2))
plot(mean(A_mse_tr,2))

%figure; hold on;
plot(mean(A0_mse_val,2))
plot(mean(A0_mse_tr,2))

%%
