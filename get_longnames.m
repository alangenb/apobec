function n = get_longnames(g)
% n = get_longnames(g)
%
% input:    g = list of gene names
% output:   n = longnames from HUGO

if isnumeric(g), error('Genenames should be strings'); end

N = load_struct('/cga/tcga-gsc/home/lawrence/db/longnames/longnames.txt');
n = mapacross(upper(g),upper(N.name),N.longname);

idx = find(strcmp('',n));
if ~isempty(idx)
  N2 = load_struct('/cga/tcga-gsc/home/lawrence/db/longnames/hugo.download.20121005.txt');
  n(idx) = mapacross(upper(g(idx)),upper(N2.ApprovedSymbol),N2.ApprovedName);
  idx = find(strcmp('',n));
  if ~isempty(idx)
    for k=1:length(idx), i=idx(k);      % now search aliases
      q1 = grepi(['(^' g{i} '$)|(, ' g{i} '$)|(^' g{i} ',)|(, ' g{i} ',)'],N2.Synonyms,1);
      if ~isempty(q1)
        n{i} = N2.ApprovedName{q1(1)};
end,end,end,end

return

%%% old method

H = load_struct('/xchip/tcga/gbm/analysis/lawrence/db/hugo.txt');
A = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');
H.gene = apply_aliases(H.ApprovedSymbol,A);
g = apply_aliases(g,A);

idx = listmap(g,H.gene);
i1 = find(~isnan(idx));
i2 = idx(i1);

n = cell(length(g),1);
n(i1) = H.ApprovedName(i2);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  make lookup table for faster processing

H = load_struct('/xchip/tcga/gbm/analysis/lawrence/db/hugo.txt');
H = reorder_struct(H,strcmp('Approved',H.Status));
N = [];
N.name = H.ApprovedSymbol;
N.longname = H.ApprovedName;
save_struct(N,'/xchip/tcga/gbm/analysis/lawrence/db/longnames.txt');
