%% Speed prediction...
% notes: 9-27-16: 5ms filter on emg doesn't seem to work...


load('./mat-files/DataNajja925_5ms');
%load('./mat-files/DataNajja925');

%% get condition idxs
c1 = Data.cond==1;
c2 = Data.cond==2;
c3 = Data.cond==3;
c4 = Data.cond==4;

%% choose data tensors
c = c1;
X = [];
X{1} = Data.Ys(:,:,c);
Y = Data.B(:,:,c);
[N,T,C] = size(Data.Ys);
P = size(Data.B,1);

%% normalize data
X{1} = bsxfun(@rdivide, X{1}, max(X{1}(:,:),[],2));
Y = bsxfun(@rdivide, Y, range(Y(:,:),2));

%% trial = 1 cycle...
cyctimes = [];
tend = 3500;
for rr = 1:size(Y,3);
  cyctimes(rr,:) = find(Y(1,1:tend-1,rr)<0 & Y(1,2:tend,rr)>0);
end

X_=[];
Y_=[];
for rr = 1:size(Y,3)
  for tt = 1:size(cyctimes,2)-1
    X_{1}(:,:,rr,tt) = X{1}(:,cyctimes(rr,tt):cyctimes(rr,tt)+450,rr);
    Y_(:,:,rr,tt) = Y(:,cyctimes(rr,tt):cyctimes(rr,tt)+450,rr);
  end
end

%% ?
X{1} = X_{1}(:,:,:);
Y = Y_(:,:,:);
%% ?

%% denoise
%mlr = [4 5 10];
mlr = 0.905;

X{2} = tensorDenoiseSVD(X_{1}, mlr);
disp(mlrank(X{2}));

X{3} = filterGauss(X{1}, 30, 2);

dat = {...
       '10ms',...
       'tensor denoised',...
       '10ms + 30ms',...
       };
   
%% regression
Yhat = [];
Xreg = [];
Yreg = Y(:,:);
for ll = 1:length(X)
  Xreg{ll} = X{ll}(:,:);
  Xreg{ll} = cat(1, Xreg{ll}, ones(1,size(Xreg{ll},2)));
  Yhat{ll} = Yreg*pinv(Xreg{ll})*Xreg{ll};
end

%% reshape
T = size(X{1},2);
for ll = 1:length(Xreg)
    Xreg{ll} = reshape(Xreg{ll}, [N+1, T, size(X{1},3)]);
    Yhat{ll} = reshape(Yhat{ll}, [P, T, size(X{1},3)]);
end

%% tensorcompare plot
rr = 1;
Z = [];
Z{1} = Y;
Z{2} = Yhat{1};
Z{3} = Yhat{2};
Z{4} = Yhat{3};
tensorComparePlot(Z,rr,3,2)
legend('data',dat{:})

%% tensorcompare plot
rr = 1;
Z = [];
Z{1} = X{1};
Z{2} = X{2};
Z{3} = X{3};
%Z{4} = Yhat{3};
tensorComparePlot(Z,rr,3,2)
%legend('data',dat{:})
legend('data','td','smooth')

%%
figure; hold all
tt = 5000;

b = 5;
plot(Yreg(b,1:tt)');
plot(Yhat{1}(b,1:tt)');
plot(Yhat{2}(b,1:tt)');
plot(Yhat{3}(b,1:tt)');
legend('data',dat{:});

%% err
err = [];
for ll = 1:length(X)
  err{ll} = sum((Yreg(:,:)-Yhat{ll}(:,:)).^2, 2)./size(Yreg(:,:),2);
end

figure
plot([err{:}]);
legend(dat{:});
ylabel('mse');
xlabel('behavioral variable: hp, vp, hv, vv, av');
%% plot
% tensorSubplot(X{1});
% tensorSubplot(X{2});
% tensorSubplot(X{3});

%% 




