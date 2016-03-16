function DataMatched = tensorDenoiseMakeDataMatched(Data, trialCount)
% to do:
% get all Y* fields.

  DataMatched = Data;
  DataMatched.Y = nan(size(Data.Y));
  DataMatched.Ys = nan(size(Data.Ys));

  for nn = 1:size(Data.Y,1)
    for cc = 1:size(Data.Y,3)
      DataMatched.Y(:,:,:,1:trialCount(nn,cc))   = Data.Y(:,:,:,1:trialCount(nn,cc));
      DataMatched.Ys(:,:,:,1:trialCount(nn,cc))  = Data.Ys(:,:,:,1:trialCount(nn,cc));
    end
  end

  DataMatched.Ysm = mean(DataMatched.Ys,4,'omitnan');
  
end



