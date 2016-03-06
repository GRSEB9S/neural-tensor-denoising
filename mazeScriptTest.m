%% RCpre
% preprocesses RC struct
clear
load('RC,2009-09-18,1-2,good-ss');
load('Nstructs_REAL');

%% 
clearvars -except R SU N*

%% trim
R = R([R.trialType]>0 & ...
      [R.unhittable]==0 & ...
      [R.novelMaze]==0 & ...
      [R.success]==1);

   pre = 400; % default: 400
   post = 100; % default: 800
   long = ([R.trialEndsTime] - [R.offlineMoveOnsetTime]) >= pre & ...
          ([R.actualDelay]) >= pre + 50; % need a buffer, i think.
           %([R.offlineMoveOnsetTime] - [R.actualFlyAppears]) >= post; % is actualFlyAppears the right field to use?
R = R(long);

%% size params
   n = size(R(1).unit,2);
   c = length(unique([R.primaryCondNum]));
   bw = 1; % bin width // try with 10
   sd = 3; % sd for guassian smoothing
   T = 1:bw:length(-pre:post);
   t = length(T);

%% spikes
   condpertrial = hist([R.primaryCondNum], unique([R.primaryCondNum]));
   maxtr = max(condpertrial);
   A = nan(n,t,maxtr,c);
   As= nan(n,t,maxtr,c);
   trxc = zeros(maxtr,c); % keep track of trials added to array
   for tr = 1:length(R)
      at = R(tr).offlineMoveOnsetTime; % alignment time
      cc = R(tr).primaryCondNum;
      thistr = find(trxc(:,cc)==0,1);
      for nn = 1:n
         % bin
         spt = floor(R(tr).unit(nn).spikeTimes)+1;
         spt = spt - (at-pre);
         spks = histc(spt,T);
         %spks = 1000/bw*histc(spt,T); % // dont use if using tensorADMM
         A(nn,:,thistr,cc) = spks;
         
         % smooth spikes
         smspks = 1000*FilterSpikes(sd,histc(spt,1:length(-pre:post)));
         As(nn,:,thistr,cc) = smspks(1:bw:length(-pre:post));
      end
      trxc(thistr,cc) = 1;
   end
   % fix later to use filterGauss
   
%% permute
   A = permute(A, [1 2 4 3]);
   As= permute(As,[1 2 4 3]);
   RC.A = A;
   RC.As = As;
   
%% load above (SKIP//)
   load('RC_10bin_10win');
   A = RC.A;
   As = RC.As;
%% 
   
%% best
   [~,is] = sort(nanmean(A(:,:),2),'descend');
   A(:,:) = A(is,:);
   As(:,:) = As(is,:);
   [n,t,c,tr] = size(A);
  
%% hosvd -- for shenoy talk
   clearvars svar rnk
   [~,S] = hosvd(nanmean(As,4));
   for ss = 1:length(S)
      sss = diag(S{ss}).^2;
      svar{ss} = cumsum(sss)/sum(sss);
      rnk(ss) = find(svar{ss} > .8 ,1); % set to 0.94 for 90% vaf
      %rnk2(ss) = find(svar{ss} > .9 , 1);
   end
   H = hosvd(nanmean(As,4),rnk);
   %% residual stuff... probably ignore
   resA = nanmean(As,4)-double(H);
   resA15 = nanmean(As,4) - nanmean(As15,4);
   resA20 = nanmean(As,4) - nanmean(As20,4);   
   for nn = 1:n
   for cc = 1:c
      xc1(nn,:,cc) = xcorr(resA(nn,:,cc),20,'coeff');
      xc2(nn,:,cc) = xcorr(resA15(nn,:,cc),20,'coeff');
      xc3(nn,:,cc) = xcorr(resA20(nn,:,cc),20,'coeff');
   end
   end
   xc1 = nanmean(nanmean(xc1,3),1);
   xc2 = nanmean(nanmean(xc2,3),1);
   xc3 = nanmean(nanmean(xc3,3),1);   
   
%%
   sortVal = 10;
   
   figA(H,1,'dots',41,'sort',sortVal)
   figA(H,2,'dots',41,'sort',sortVal)
   figA(H,3,'dots',41,'sort',sortVal)
   figA(H,4,'dots',41,'sort',sortVal)
   %%
   figA(H1,1,'dots',41,'sort',sortVal)
   figA(H1,2,'dots',41,'sort',sortVal)
   figA(H1,3,'dots',41,'sort',sortVal)
   figA(H1,4,'dots',41,'sort',sortVal)
   %%
   figA(H2,1,'dots',41,'sort',sortVal)
   figA(H2,2,'dots',41,'sort',sortVal)
   figA(H2,3,'dots',41,'sort',sortVal)
   figA(H2,4,'dots',41,'sort',sortVal)
   %%
   figA(nanmean(As,4),1,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),2,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),3,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),4,'dots',41,'sort',sortVal)
   %%
   figA(nanmean(As,4),1,'dots',41,'sort',sortVal)
   
   figA(H,2,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),2,'dots',41,'sort',sortVal)
   figA(H,3,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),3,'dots',41,'sort',sortVal)
   figA(H,4,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),4,'dots',41,'sort',sortVal)
%%
   recerr1 = 1-sum(vec((double(H)  - nanmean(As,4)).^2))/sum(vec(nanmean(As,4)).^2)
   recerr2 = 1-sum(vec((double(H1) - nanmean(As,4)).^2))/sum(vec(nanmean(As,4)).^2)   
   recerr3 = 1-sum(vec((double(H2) - nanmean(As,4)).^2))/sum(vec(nanmean(As,4)).^2)
%%
   recerr4 = 1-sum(vec((nanmean(As25,4) - nanmean(As,4)).^2))/sum(vec(nanmean(As,4)).^2)
%%
   figA(nanmean(As,4),1,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),2,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),3,'dots',41,'sort',sortVal)
   figA(nanmean(As,4),4,'dots',41,'sort',sortVal)
%% tensor denoising
   nd = 100;
   td = t;
   cd = 100;

%%
%% ------ FOR VISUALIZATION
%% ------ do tensorADMM, then hosvd to remove very low SVs... tensorADMM might always make things noisy because its never removing noise, just making it small
T1 = tensorADMM(A(1:nd,1:td,1:cd,:),'ef',1,'bf',.1,'gf',1,'s',1,'noiseTerm','poiss','maxIter',50);

% for M1 data: ef 1 bf 0.1 gf 1 s 1 works really well (I think)
%%
T1{1} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',1,'bf',1','s',1,'sig',0.01,'noiseTerm','gauss','maxIter',10);
T1{2} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',1,'bf',1','s',1,'sig',0.1,'noiseTerm','gauss','maxIter',10);
T1{3} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',1,'bf',1','s',1,'sig',1,'noiseTerm','gauss','maxIter',10);

T2{1} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',10,'bf',1','s',1,'sig',0.1,'noiseTerm','gauss','maxIter',10);
T2{2} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',1,'bf',1','s',1,'sig',0.1,'noiseTerm','gauss','maxIter',10);
T2{3} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',.1,'bf',1','s',1,'sig',0.1,'noiseTerm','gauss','maxIter',10);

T3{1} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',1,'bf',10','s',1,'sig',0.1,'noiseTerm','gauss','maxIter',10);
T3{2} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',1,'bf',1','s',1,'sig',0.1,'noiseTerm','gauss','maxIter',10);
T3{3} = tensorADMM(As(1:nd,1:td,1:cd,:),'ef',1,'bf',.1','s',1,'sig',0.1,'noiseTerm','gauss','maxIter',10);












