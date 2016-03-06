%% tensor ADMM demo
%% build underlying rates
   clear; close all; clc
   n =  40; %neurons
   nd = 10; %neural dimensions
   c =  40;
   cd = 10;
   t = 80;
   td = 30;

   tr = 50; %trials
   
   %[~,~,~,xnat,~] = RDLTI(nd,t,c,1,n,cd);
   
   
   
   % nd x td x cd tensor
   xnat = randn(nd,td,cd); % core tensor.
   W{1} = orth(randn(n,nd)); % left SVs for x_(1)
   % patterns, left SVs for x_(2)
      tps = exp((-0.05*rand(td,1) + 0.1*randn(td,1)*1i)*(1:t));
      tps = bsxfun(@times, exp(2*pi*1i*randn(td,1)), tps);
      W{2} = real(tps)';
   W{3} = orth(randn(c,cd)); % left SVs for x_(3)
   xnat = 2.*double(ttm(tensor(xnat),W)); % x_natural, rates before passing through exp()

   %
   offset = round(t/8);
   stp = zeros(size(xnat));
   stp(:,offset:end-offset,:) = 1;
   stp = filterGauss(stp,1,2);
   xnat = xnat.*stp + bsxfun(@times, 0.1*rand(n,1), stp);
   xnat = bsxfun(@plus, randn(n,1), xnat);
   
   s = 1;
   f = @(x,s) s*log(1 + exp(x/s));
   x = 0.5.*f(xnat,s);
   
%% poisson observations
   Y = zeros(n,t,c,tr);
   for rr = 1:tr
      Y(:,:,:,rr) = random('poiss',x);
   end
   disp(mean(Y(:)))
%% -- add smoothing Ys, test with gaussian noise
   Ys = filterGauss(Y,2,2);

%%
T1 = tensorADMM(Y,'ef',1,'bf',.1,'gf',1,'s',1,'noiseTerm','poiss','maxIter',150);
%%
[T B] = tensorADMM2(Y,'ef',1,'bf',.3,'gf',1e2,'s',1,'noiseTerm','poiss','maxIter',100);

%%
figure(10);
for i = 1:30;
    for ii = 1:5
    subplot(5,1,ii)
    plot(T3c(ii,:,i));
    hold on;plot(log(1+exp(repmat(Bc(ii),1,t,1)+T3c(ii,:,i))),'r--');
    hold off;
    end
    pause
end
%%
plotA(xnat)
plotA(T1);
plotA(mean(Y,4));

%%
plotAgrid(T1(:,:,1:10),3,3);
plotAgrid(xnat(:,:,1:10),3,3)
Ym = mean(Y,4);
plotAgrid(Ym(:,:,1:10),3,3);









