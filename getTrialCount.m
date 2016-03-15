function trialCount = getTrialCount(Yin)
% trialCount = getTrialCount( Yin )
% takes a size (neuron,time,condition,trial) tensor and counts how many
% non-nan entries are in each (i,1,j,:) vector. 

[n,t,c,r] = size(Yin);
trialCount = zeros(n,c);
  for nn = 1:n
    for cc = 1:c
      trialCount(nn,cc) = sum(~isnan(Yin(nn,1,cc,:)));
    end
  end
end
