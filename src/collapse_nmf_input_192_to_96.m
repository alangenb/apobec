function X = collapse_nmf_input_192_to_96(X)

demand_field(X,{'pat','chan'});
demand_field(X.pat,'nchan');

if size(X.pat.nchan,2)==96, fprintf('Already collapsed!\n'); return; end
if size(X.pat.nchan,2)~=192, error('Input is not 192-wide!\n'); end

X.chan = parsein(X.chan,'name','^(.) \((.)->(.)\) (.)$',{'l','f','t','r'});
X.chan.l = listmap(X.chan.l,{'A','C','G','T'});
X.chan.f = listmap(X.chan.f,{'A','C','G','T'});
X.chan.t = listmap(X.chan.t,{'A','C','G','T'});
X.chan.r = listmap(X.chan.r,{'A','C','G','T'});

keep = true(slength(X.chan),1);
rc = [4 3 2 1];
for f=1:2, for t=1:4, for l=1:4, for r=1:4
  if f==t, continue; end
  idx1 = find(X.chan.l==l & X.chan.f==f & X.chan.t==t & X.chan.r==r);
  idx2 = find(X.chan.l==rc(r) & X.chan.f==rc(f) & X.chan.t==rc(t) & X.chan.r==rc(l));
  if length(idx1)~=1 || length(idx2)~=1, error('problem with X.chan'); end
  X.pat.nchan(:,idx1) = X.pat.nchan(:,idx1) + X.pat.nchan(:,idx2);
  keep(idx2) = false;
end,end,end,end

X.chan = reorder_struct(X.chan,keep);
X.pat.nchan = X.pat.nchan(:,keep);

