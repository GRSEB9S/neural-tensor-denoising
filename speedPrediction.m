%% Speed prediction...

load('./mat-files/DataNajja');

%% get condition idxs
c1 = T1.aligned.conditionID==1;
c2 = T1.aligned.conditionID==2;
c3 = T1.aligned.conditionID==3;
c4 = T1.aligned.conditionID==4;

%% choose data tensors
c = c1;
X = [];
X{1} = Data.Ys(:,:,c);
Y = Data.B(:,:,c);

%% normalize data
X{1} = bsxfun(@rdivide, X{1}, max(X{1}(:,:),[],2));
Y = bsxfun(@rdivide, Y, range(Y(:,:),2));

%% denoise
mlr = [5 10 20];
%mlr = 0.90;
X{2} = tensorDenoiseSVD(X{1}, mlr);
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
  Yhat{ll} = Yreg/Xreg{ll}*Xreg{ll};
end


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
  err{ll} = sum((Yreg-Yhat{ll}).^2, 2)./size(Yreg,2);
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




