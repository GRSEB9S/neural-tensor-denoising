function [rnks,x,y] = paretoPoints(rnks, err, coreMeasure)
%% pareto optimal
% get pareto optimal points in the set of core tensor sizes, for a given
% choice of model complexity measure.

  x = coreMeasure(rnks);
  y = err;

  d = sqrt(x.^2+y.^2);
  [d,idx] = sort(d);
  rnks = rnks(idx,:);
  x = x(idx);
  y = y(idx);

  p = 1;
  while sum(p) <= length(x) % originally it was ~=. now it is <=
      idx = (x >= x(p) & y > y(p)) | (x > x(p) & y >= y(p));
      if any(idx)
          rnks = rnks(~idx,:); x = x(~idx); y = y(~idx); d = d(~idx);
      end
      p = p+1-sum(idx(1:p));
  end
  [~,idx] = sortrows([coreMeasure(rnks) -y]);
  rnks = rnks(idx,:); x = x(idx); y = y(idx);
end
