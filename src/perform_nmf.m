function X = perform_nmf(X,k,randseed)

if nargin<2, error('usage: perform_nmf(X,k)  where X is from make_input_for_nmf'); end
if ~exist('randseed','var'), randseed=1234; end

rand('twister',randseed);
[X.pat.nmf,X.chan.nmf]=nmf_wrapper(X.pat.nchan,k);

X.chan.nmf = X.chan.nmf';
X.pat.nmf_norm = bsxfun(@rdivide,X.pat.nmf,nansum(X.pat.nmf,2));

% compute # mutations assigned to each process
X.pat.nmf_nmut = nan(slength(X.pat),k);
for i=1:k
  X.pat.nmf_nmut(:,i) = nansum(X.pat.nmf(:,i)*X.chan.nmf(:,i)',2);
end




