%% sfn2014_4 -- trial reduction using simulated data
%% load data
load('sfn_10bw_10sd.mat');
%% prepare arrays
A  = D.L.A;
T  = D.L.T;
T0 = D.L.T0;
bw = D.bw;
[n,t,c,r] = size(T);
bf = bw/1000;

%%
s = 0.1;
f = @(x,s) s*log(1+exp(x/s));
%% create signals
A_tens_true = double(tensorOpt(nanmean(A,4), 50, 1));
A_tens_true = A_tens_true.*(A_tens_true>0);
A_tens_true = bf*A_tens_true;

R = 100;
Y = zeros(n,t,c,R);
for rr = 1:R
   Y(:,:,:,rr) = poissrnd(A_tens_true);
end
%% params
P = 4;
A0param = linspace(2,5,P);

%% increase R
for rr = 1:R
   disp(['Tensor model. r = ' num2str(rr)]);
   A_tens{rr} = double(tensorOpt(nanmean(Y(:,:,:,1:rr),4),50,1));
end
%%
for rr = 1:R
   disp(['Filter model. r = ' num2str(rr)]);
   for pp = 1:P
      A0_tens{pp,rr} = nanmean(filterGauss(Y(:,:,:,1:rr),A0param(pp),2),4);
   end
end

%% compute errors
for rr = 1:R
   disp(['Tensor error. r = ' num2str(rr)]);
   
   %A_ll(rr)  =  sfn_ll(Y(:,:,:,1:rr), f(A_tens{rr},s));
   A_mse(rr) = sfn_mse(A_tens_true, A_tens{rr});

end
%%
for rr = 1:R
   disp(['Filter error. r = ' num2str(rr)]);
   
   for pp = 1:P
      A0_mse(pp,rr) = sfn_mse(A_tens_true, A0_tens{pp,rr});
   end

end

%%



