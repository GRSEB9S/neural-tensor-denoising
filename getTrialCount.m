
function trialOut = getTrialCount(Yin)
  [n,t,c,r] = size(Yin);
  for nn = 1:n
    for cc = 1:c
      trialOut(nn,cc) = sum(~isnan(Yin(nn,1,cc,:)));
    end
  end
end