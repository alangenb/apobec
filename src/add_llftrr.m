function m = add_llftrr(m,strand_collapse,path)
% given a mutation list, adds the simple ll,l,f,t,r,rr (1-4) annotations
% if strand_collapse==1, flips G->C and T->A

if ~exist('strand_collapse','var'), strand_collapse = true; end

if isfield(m,'newbase') && iscellstr(m.newbase)
  m.t = listmap(m.newbase,{'A','C','G','T'});
elseif isfield(m,'newbase_idx')
  m.t = m.newbase_idx;
elseif isfield(m,'alt_idx')
  m.t = m.alt_idx;
elseif isfield(m,'alt') && isnumeric(m.alt) && mean(m.alt>=1 & m.alt<=4)>=0.8
  m.t = m.alt;
elseif isfield(m,'alt') && iscellstr(m.alt)
  m.t = listmap(m.alt,{'A','C','G','T'});
else
  m.newbase = find_newbase(m);
  m.t = listmap(m.newbase,{'A','C','G','T'});
end

if ~isfield(m,'context1025')
  demand_fields(m,{'chr','pos'});
  m.context1025 = get_context1025(m.chr,m.pos,path);
end

d = load_categs('context1025');
d = parsein(d,'name','^(.) in (.)(.)_(.)(.)$',{'f','ll','l','r','rr'});
d.f = listmap(d.f,{'A','C','G','T'});
d.l = listmap(d.l,{'A','C','G','T'});
d.r = listmap(d.r,{'A','C','G','T'});
d.ll = listmap(d.ll,{'A','C','G','T'});
d.rr = listmap(d.rr,{'A','C','G','T'});
m = assigninto(m,d,m.context1025);
m = mf2a(m,'t','f');

% strand-collapse if requested
if strand_collapse
  fprintf('Collapsing strands.\n');
  gt = (m.f>2);
  m.f(gt) = 5-m.f(gt);
  m.t(gt) = 5-m.t(gt);
  tmp = 5-m.l(gt);
  m.l(gt) = 5-m.r(gt);
  m.r(gt) = tmp;
  tmp = 5-m.ll(gt);
  m.ll(gt) = 5-m.rr(gt);
  m.rr(gt) = tmp;
else
  fprintf('Not collapsing strands.\n');
end

% identify non-SNPs if possible, NaN them out
bad = (isnan(m.f)|isnan(m.t));
if isfield(m,'classification') bad(~strcmp('SNP',m.classification))=1; end
if isfield(m,'type') bad(grepmi('ins|del|frame',m.type))=1; end
m.f(bad)=nan;
m.l(bad)=nan;
m.r(bad)=nan;
m.t(bad)=nan;
m.ll(bad)=nan;
m.rr(bad)=nan;







