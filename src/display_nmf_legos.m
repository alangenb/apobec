function display_nmf_legos(X,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'display','counts');
P=impose_default_value(P,'fontsize',16);

if isfield(X,'chan')
  data = X.chan;
elseif isfield(X,'nmf') && isfield(X,'catnames')
  data = X;
end


q={}; for i=1:size(data.nmf,2),q{i} = data.nmf(:,i);end
P.catnames=data.catnames;

if isfield(X,'factor') && isfield(X.factor,'name')
  titles = X.factor.name;
else
  titles=num2cellstr(1:length(q));
end

legos(q,titles,P);ff
