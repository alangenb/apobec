function [rho pval] = corrr(x,y)

x=double(as_column(x));
y=double(as_column(y));

idx = find(~isnan(x)&~isnan(y));
if ~isempty(idx)
  [rho pval] = corr([x(idx) y(idx)]);
  rho = rho(1,2);
  pval = pval(1,2);
else
  rho = nan;
  pval = nan;
end

