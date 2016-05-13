%% load, use tensorDenoiseMakeData()

%% analyze
tensorSubplot(mean(DataA.Ys(:,:,:,1:2),4))
tensorSubplot((DataA.Ysm))
tensorSubplot(tensorDenoiseSVD(mean(DataA.Ys(:,:,:,1:2),4),0.975));


%%
















