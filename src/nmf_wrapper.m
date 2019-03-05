function [w h] = nmf_wrapper(varargin)

v = varargin{1};
k = varargin{2};

% REMOVE all-NaN rows/columns
xidx = (~all(v==0 | isnan(v),1));
yidx = (~all(v==0 | isnan(v),2));
v2 = v(yidx,xidx);

% NMF 
[w1 h1] = nmf(v2, k, varargin{3:end});

% SORT the <k> components so they are always in the same order no matter what random seed was used
%  (1) choose the top N highest-variance elements (in the longer dimension)
%  (2) sort the k components by these elements
if size(v2,1)>size(v2,2)

  ve = sum(bsxfun(@times,w1,var(v2,0,2)),1);
  [tmp ord] = sort(ve,'descend');

%  ncomp = 20;
%  q = var(v2,0,2);
%  [tmp ord] = sort(q,'descend');
%  w2 = w1(ord(1:ncomp),:);
%  for i=1:size(w2,1)
%    mx=max(w2(i,:));
%    w2(i,w2(i,:)<mx)=0;
%  end
%  [tmp ord] = sortrows(w2'); ord=flip(ord);

else % sort based on cols of v
  ve = sum(bsxfun(@times,h1,var(v2,0,1)),2);
  [tmp ord] = sort(ve,'descend');
end

w1 = w1(:,ord);
h1 = h1(ord,:);

% REPLACE NaNs
xmap = nan(size(v,2),1);
ymap = nan(size(v,1),1);
xmap(xidx)=1:size(h1,2);
ymap(yidx)=1:size(w1,1);
w = nansub(w1,ymap);
h = nansub(h1',xmap)';

