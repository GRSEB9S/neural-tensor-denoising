function Yout = resampleTrials( Yin, repl, seed )
  
  if nargin < 2
    repl = 1;
  end
  if nargin < 3
    seed = [];
  end
    
  if ~isempty(seed)
    rng(seed,'twister');
  end
  
  trialCount = getTrialCount(Yin);
  [n,t,c,r] = size(Yin);
  for nn = 1:n;
    for cc = 1:c
      % replacement option 1: with replacement
      if repl
        howmany = r;
        trialInds = randi(trialCount(nn,cc), 1, howmany);
      % replacement option 2: without replacement
      else
        howmany = trialCount(nn,cc);
        trialInds = randperm(trialCount(nn,cc));
      end
      Yin(nn,:,cc,1:howmany) = Yin(nn,:,cc,trialInds);
    end
  end
  Yout = Yin;
end