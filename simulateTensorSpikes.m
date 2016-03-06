function [ Yout, Xsim, Bsim, s, f ] = simulateTensorSpikes
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

   n =  50; %neurons
   nd = 10; %neural dimensions / bases
   c =  30; %conditions (e.g. stimuli)
   cd = 6; %condition dimensions / bases
   t = 300;  %time bins
   td = 6; %temporal dimensions / bases
   
   Xsim = 2*trnd(10,nd,td,cd); % start with core tensor.
   % bases for neurons / left SVs for x_(1)
      W{1} = orth(randn(n,nd)); % orth / no orth?
   % temporal patterns / left SVs for x_(2)
      tbuffer = 20;
      t = t+tbuffer;
      tps = exp((2*pi/t*(1.5*randn(td,1) + 0))*1i*(1:t));
      tps = bsxfun(@times, exp(2*pi*1i*rand(td,1)), tps);
      W{2} = qr(real(tps))';
      W{2} = W{2}(tbuffer+1:end,:);
      t = t-tbuffer;
   % bases for conditions / left SVs for x_(3)
      W{3} = orth(randn(c,cd)); % orth / no orth?
   
   Xsim = 1*double(ttm(tensor(Xsim),W)); % build up tensor using ttm (tensor times matrices).
   Bsim = 1*randn(n,1) - 2;
   
   s = 1;
   f = @(x,b,s) s*log(1 + exp(bsxfun(@plus, x, b)/s)); % firing rate function

  %% spikes
  r = 40; %trials
  Ysim = zeros(n,t,c,r);
  for rr = 1:r
    Ysim(:,:,:,rr) = poissrnd(0.05*f(Xsim,Bsim,s));
  end
  disp(mean(vec(mean(Ysim,4))))

  %% lopsided
  c1 = round(0.1*c);
  r1 = round(0.1*r);

  Ysim2 = Ysim;
  Ysim2(:,1:c1,:,r1+1:end) = nan;
  
  %% missing data
  Ysim3 = Ysim;
  c1 = round(0.1*c);
  for nn = 1:n
    cp = randperm(c);
    Ysim3(nn,:,cp(c1+1:end)) = nan;
  end

  %%
  Yout.Y = Ysim;
  Yout.Yl = Ysim2;   
  Yout.Ym = Ysim3; 
  
end






















