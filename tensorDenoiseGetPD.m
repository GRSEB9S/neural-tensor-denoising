function pd = tensorDenoiseGetPD(Y, targ)

  [n,t,c] = size(Y);
  pd = zeros(n,t);
  for nn = 1:n
    for tt = 1:t
      tc = squeeze(Y(nn,tt,:))';
      pvec = tc/targ;
      pd(nn,tt) = atan2(pvec(1), pvec(2))*180/pi;
    end
  end

end