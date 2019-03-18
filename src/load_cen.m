function cen = load_cen(build)
if ~exist('build','var'), build='hg18'; end
B=load_cytobands(build);
c=grep('cen',B.stain,1);
cen = zeros(24,2);
for i=1:24
  idx = intersect(c, find(B.chr==i));
  cen(i,1) = min(B.start(idx));
  cen(i,2) = max(B.end(idx));
end
