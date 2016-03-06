%% tensor ADMM demo, for sfn abstract 
%  never actually used for sfn abstract
%% build underlying rates
   clear; close all; clc
   n =  50; %neurons
   nd = 10; %neural dimensions
   c =  50;
   cd = 10;
   t = 50;
   td = 5;

   tr = 100; %trials

   % nd x td x cd tensor
   xnat = randn(nd,td,cd); % core tensor.
   W{1} = orth(randn(n,nd)); % left SVs for x_(1)
   % patterns, left SVs for x_(2)
      for tt = 1:td
         W{2}(:,tt) = sin(pi*tt/t*(0:t-1) + pi*rand);
      end
      
      
      
      %tps = randn(t+100,1000);
      %tps = filter(gausswin(t/2),1,tps);
      %[tps,~,~] = svd(tps(end-t+1:end,:));
      %W{2} = tps(:,1:td);
   W{3} = orth(randn(c,cd)); % left SVs for x_(3)
   xnat = 1*double(ttm(tensor(xnat),W)); % x_natural, rates before passing through exp()

   f = @(x) log(1 + exp(x));
   %x = exp(xnat); % nonnegative rate exp(xnat)
   x = .1*f(xnat);
   
%% poisson observations
   Y = zeros(n,t,c,tr);
   for rr = 1:tr
      Y(:,:,:,rr) = random('poiss',x);
   end
   Ym = mean(Y,4);
   mean(Y(:)) % avg spikes per bin, should be ~ 0.07 to match m1
   
   
   
   Ys = Y;
   for nn = 1:n
   for cc = 1:c
   for rr = 1:tr
      Ys(nn,:,cc,rr) = FilterSpikes(1.5,Y(nn,:,cc,rr));
   end
   end
   end
   
%%
tr_1 = 20;
Xadmm1 = tensorADMM(Y(:,:,:,1:tr_1),xnat,1,.3,.3);
Xadmm2 = tensorADMM(Y(:,:,:,1:tr_1),xnat,1,.3,.03);
Xadmm3 = tensorADMM(Y(:,:,:,1:tr_1),xnat,1,.1,.03);
Xadmm4 = tensorADMM(Y(:,:,:,1:tr_1),xnat,1,.1,.01);
%Xadmm5 = tensorADMM(Y(:,:,:,1:tr_1),xnat,1,10,1);

%%

Ym2 = mean(Y(:,:,:,1:tr_1),4);

Ys2 = Ym2;
for nn = 1:n
for cc = 1:c
   Ys2(nn,:,cc) = FilterSpikes(2,Ym2(nn,:,cc));
end
end


%%
xe_ten = 100*x - 100*f(Xadmm1);
xe_ten = xe_ten.^2;
xe_ten = mean(xe_ten(:))

%%
ye_mn = 100*x - 100*mean(Y(:,:,:,1:tr_1),4);
ye_mn = ye_mn.^2;
ye_mn = mean(ye_mn(:))

%%
ye_mn2 = 100*x - 100*Ym2;
ye_mn2 = ye_mn2.^2;
ye_mn2 = mean(ye_mn2(:))

%%
ye_mn2s = 100*x - 100*Ys2;
ye_mn2s = ye_mn2s.^2;
ye_mn2s = mean(ye_mn2s(:))


%%
for rr = 1:tr
   ym = 100*mean(Y(:,:,:,1:rr),4);
   ym = 100*x - ym;
   ym = ym.^2;
   ye_mn2(rr) = mean(ym(:));
end
%%
figure; hold all;

plot(ye_mn);
plot([1 tr],[xe_ten xe_ten]);
plot([1 tr],[ye_mn ye_mn]);



%%
nc = 4;
%%
plotAgrid(x(:,:,1:nc),2,2)

plotAgrid(f(Xadmm1(:,:,1:nc)),2,2)
plotAgrid(f(Xadmm2(:,:,1:nc)),2,2)
plotAgrid(f(Xadmm3(:,:,1:nc)),2,2)
plotAgrid(f(Xadmm4(:,:,1:nc)),2,2)

%%
plotAgrid2(f(Xadmm1(:,:,1:nc)),x(:,:,1:nc),2,2);
plotAgrid2(f(Xadmm2(:,:,1:nc)),x(:,:,1:nc),2,2);
plotAgrid2(f(Xadmm3(:,:,1:nc)),x(:,:,1:nc),2,2);
plotAgrid2(f(Xadmm4(:,:,1:nc)),x(:,:,1:nc),2,2);

%%
plotAgrid2(Ys2(:,:,1:nc),x(:,:,1:nc),2,2);
plotAgrid2(Ym2(:,:,1:nc),x(:,:,1:nc),2,2);



%%














