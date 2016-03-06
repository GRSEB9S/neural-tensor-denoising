%% for AI poster
%% load data
load('sfn_10bw_10sd.mat');
%% prepare arrays
A  = D.L.A;
T  = D.L.T;
T0 = D.L.T0;
bw = D.bw;
[n,t,c,r] = size(A);
bf = bw/1000;
Ablah = squeeze(A(:,1,:,:));
Ablah = permute(Ablah,[3 2 1]);
nnnn = sum(~isnan(Ablah(:,:)));
%% fig properties
% move is at t = 500 for 1300 times
% move is at t = 51 for 130 times
figap = {'dots',51,...
          'ps',[100 40],...
          'ybar',20,...
          'xbar',20,...
          'csort',0,...
          'cmap','lara',...
          };
       
%%
%time = 30:100;
time = 1:t;
%%
f_ = @(x) x.*(x>0);
%% panel 1
% start with noisy histograms (low trials) -> smoothing 
r = 5;
Aa = nanmean(A(:,:,:,1:2),4);
Ab = nanmean(A,4);
Ac = nanmean(filterGauss(Ab,2.5,2),4);

Ad = filterGauss(double(tensorOpt(Aa, 35, 1)),1.5,2);
Add = filterGauss(double(tensorOpt(Ab, 63, 1)),1.5,2);
%%
% 97, 23, 37 (87, 10)...101?
       
nn = 97;

figA(Aa(:,time,:),nn,figap{:}); h(1) = gca;
figA(Ab(:,time,:),nn,figap{:}); h(2) = gca;
figA(Ac(:,time,:),nn,figap{:}); h(3) = gca;
figA(f_(Ad(:,time,:)),nn,figap{:}); h(4) = gca;

y = max([h(:).YLim]);
for xx = 1:4
   h(xx).YLim = [0 y];
end

%% tensor denoising cascade
figA(Aa(:,time,:),23,figap{:});
figA(Aa(:,time,:),37,figap{:});
figA(Aa(:,time,:),87,figap{:});
figA(Aa(:,time,:),10,figap{:});


%% panel 4 -- many psths
%% check
nn = 97;
figA(Ab(:,time,:),nn,figap{:});h(1) = gca;
figA(Ac(:,time,:),nn,figap{:});h(2) = gca;
figA(f_(Add(:,time,:)),nn,figap{:});h(3) = gca;
y = max([h(1:3).YLim]);
for xx = 1:3
   h(xx).YLim = [0 y];
end
%%
nn = 23;
figA(Ab(:,time,:),nn,figap{:});h(1) = gca;
figA(Ac(:,time,:),nn,figap{:});h(2) = gca;
figA(f_(Add(:,time,:)),nn,figap{:});h(3) = gca;
y = max([h(1:3).YLim]);
for xx = 1:3
   h(xx).YLim = [0 y];
end
%%
nn = 37;
figA(Ab(:,time,:),nn,figap{:});h(1) = gca;
figA(Ac(:,time,:),nn,figap{:});h(2) = gca;
figA(f_(Add(:,time,:)),nn,figap{:});h(3) = gca;
y = max([h(1:3).YLim]);
for xx = 1:3
   h(xx).YLim = [0 y];
end
%%
nn = 87;
figA(Ab(:,time,:),nn,figap{:});h(1) = gca;
figA(Ac(:,time,:),nn,figap{:});h(2) = gca;
figA(f_(Add(:,time,:)),nn,figap{:});h(3) = gca;
y = max([h(1:3).YLim]);
for xx = 1:3
   h(xx).YLim = [0 y];
end

%% panel 5 -- validation
load('_3', 'A_mse_val', 'A_mse_tr', 'Aparam');
figure; hold on;
plot(Aparam,mean(A_mse_val/100^2,2))
plot(Aparam,mean(A_mse_tr/100^2,2))
[~,mi] = min(mean(A_mse_val,2));
rnk = Aparam(mi);
plot(rnk, mean(A_mse_val(mi,:),2)/100^2,'o')
xlabel('Tensor rank');
ylabel('MSE');
ylim([0 0.02])
xlim([3 120])

%% fig examples
A1 = filterGauss(double(tensorOpt(Ab, Aparam(8), 1)),1.5,2);
A2 = filterGauss(double(tensorOpt(Ab, Aparam(mi), 1)),1.5,2);
A3 = filterGauss(double(tensorOpt(Ab, Aparam(end), 1)),1.5,2);
%%
nn = 97;
figA(f_(A1(:,time,:)),nn,figap{:}); h(1) = gca;
figA(f_(A2(:,time,:)),nn,figap{:}); h(2) = gca;
figA(f_(A3(:,time,:)),nn,figap{:}); h(3) = gca;
figA(Ab(:,time,:),nn,figap{:}); h(4) = gca;
y = max([h(1:4).YLim]);
for xx = 1:4
   h(xx).YLim = [0 y];
end


%% panel 6 -- trial reduction
load('_4', 'A_tens_true', 'A_tens', 'A0_tens', 'Y', 'A_mse', 'A0_mse');
figure; hold on;
xlabel('trials');
ylabel('MSE');
s = 3;
plot(A_mse)
plot(A0_mse(s,:))
xlim([1 80])
ylim([0 0.003])

tr1 = 7;
tr2 = 10;

v1 = A_mse(tr1);
tr_1 = find(A0_mse(s,:)<v1,1);
v2 = A0_mse(s,tr_1);

w1 = A_mse(tr2);
tr_2 = find(A0_mse(s,:)<w1,1);
w2 = A0_mse(s,tr_2);

plot([tr1 tr_1 tr_1],[v1 v1 0],'k');
plot([tr2 tr_2 tr_2],[w1 w1 0],'k');

%% ground truth cascade
figA(100*filterGauss(f_(A_tens_true(:,time,:)),2.5,2),23,figap{:});
figA(100*filterGauss(f_(A_tens_true(:,time,:)),2.5,2),37,figap{:});
figA(100*filterGauss(f_(A_tens_true(:,time,:)),2.5,2),87,figap{:});
figA(100*filterGauss(f_(A_tens_true(:,time,:)),2.5,2),10,figap{:});
%% example figs
% good ex: n 1
nn = 1;
figA(100*filterGauss(A_tens_true,2.5,2),nn,figap{:}); h(1) = gca;
%% Y
figA(100*filterGauss(nanmean(Y(:,:,:,1:7),4),1,2),nn,figap{:}); h(2) = gca;
figA(100*filterGauss(nanmean(Y(:,:,:,1:17),4),1,2),nn,figap{:}); h(3) = gca;
figA(100*filterGauss(nanmean(Y(:,:,:,1:100),4),1,2),nn,figap{:}); h(4) = gca;
%% A_tens
figA(100*filterGauss(f_(A_tens{7}),1.5,2),nn,figap{:}); h(5) = gca;
figA(100*filterGauss(f_(A_tens{17}),1.5,2),nn,figap{:}); h(6) = gca;
figA(100*filterGauss(f_(A_tens{100}),1.5,2),nn,figap{:}); h(7) = gca;

%% A0_tens
figA(100*A0_tens{3,7},nn,figap{:}); h(8) = gca;
figA(100*A0_tens{3,17},nn,figap{:}); h(9) = gca;
figA(100*A0_tens{3,100},nn,figap{:}); h(10) = gca;
y = max([h([1 5:10]).YLim]);
for xx = [1 5:10]
   h(xx).YLim = [0 y];
end

%%
nn = 13;
cc = 5;
Atmp = cat(3, ...
           A_tens_best(:,:,cc),...
           A0_tens_best(:,:,cc),...
           nanmean(A(:,:,cc,6:end),4),...
           nanmean(A(:,:,cc,1:5),4));
figA(Atmp(:,time,:),nn,figap{:},'cmap','lines');


%% panel 7 -- low trial reconstruction
load('_5', 'A_tens_best', 'A0_tens_best', 'A_mse_test', 'A0_mse_test');
%figure; hold on;
% % good neuron? 5, 18? 100?
nn = 103;
cc = 1:24;
time= 2:t;
figA(nanmean(A(:,time,cc,1:3),4),nn,figap{:});
figA(nanmean(A(:,time,cc,1:end),4),nn,figap{:});
Acc = nanmean(filterGauss(Ab,2,2),4); figA(Acc(:,time,cc),nn,figap{:});
figA(A_tens_best(:,time,cc),nn,figap{:});
figA(A0_tens_best(:,time,cc),nn,figap{:});

%%
figure; hold on;
% ex; n 5 c 13
%     n 68 c 10
cc = 13;
nn = 5;
plot(A_tens_best(nn,time,cc),'linewidth',2);%,nn,figap{:});
plot(A0_tens_best(nn,time,cc),'--','linewidth',2);%,nn,figap{:});
plot(nanmean(A(nn,time,cc,6:end),4));%,nn,figap{:})
plot(nanmean(A(nn,time,cc,1:5),4));%,nn,figap{:});

%%
cc = 13;
Atmp = cat(3, ...
           A_tens_best(:,:,cc),...
           A0_tens_best(:,:,cc),...
           nanmean(A(:,:,cc,6:end),4),...
           nanmean(A(:,:,cc,1:5),4));
figA(Atmp(:,time,:),nn,figap{:},'cmap','lines');
        
        
%%
bar([0 1],[A_mse_test A0_mse_test]);




