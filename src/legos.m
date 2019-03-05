function legos(varargin)
% legos(input1,input2,input3)
% legos({input1,input2,input3},titles)
% legos({input1,input2,input3},titles,P)

z = varargin{1};

if isnumeric(z) && ndims(z)==3 && size(z,1)==8 && size(z,2)==12 && size(z,3)>1
  % convert multi-page 8x12 lego matrix into cell-format
  tmp=z; z={}; for i=1:size(tmp,3),z{i,1}=tmp(:,:,i);end
end

if nargin==1 && iscell(z)
  titles = [];
  P=[];
elseif nargin==2 && iscell(z) && iscell(varargin{2})
  titles = varargin{2};
  P=[];
elseif nargin==2 && iscell(z) && isstruct(varargin{2})
  titles = repmat({''},length(varargin{1}));
  P=varargin{2};
elseif nargin==3 && iscell(z) && iscell(varargin{2}) && (ischar(varargin{3})&&strcmp('genome',varargin{3}))
  titles = varargin{2};
  P = [];
  P.display = 'rates'; % was missing until 2019-02-24
  P.coverage = 'genome';
elseif nargin==3 && iscell(z) && iscell(varargin{2}) && (ischar(varargin{3})&&strcmp('exome',varargin{3}))
  titles = varargin{2};
  P = [];
  P.display = 'rates'; % was missing until 2019-02-24
  P.coverage = 'exome';
elseif nargin==3 && iscell(z) && iscell(varargin{2}) && (ischar(varargin{3})&&strcmp('counts',varargin{3}))
  titles = varargin{2};
  P = [];
  P.display = 'counts';
elseif nargin==3 && iscell(z) && iscell(varargin{2}) && isstruct(varargin{3})
  titles = varargin{2};
  P = varargin{3};
else
  z = varargin;
  titles = [];
  P=[];
end

P = impose_default_value(P,'nrows',[]);
P = impose_default_value(P,'fontsize',12);

n = length(z);
if ~isempty(P.nrows)
  nrows=P.nrows;
  ncols = ceil(n/nrows);
else
  if n<=3
    nrows=1;
    ncols=n;
  else
    s = ceil(sqrt(n));
    nrows=s;
    ncols=s;
    while((nrows-1)*ncols>=n) nrows=nrows-1; end
  end
end

fi=1;
for y=1:nrows, for x=1:ncols
    if fi<=n
      subplot('position',[(x-1)*(1/ncols) 1-(y*(1/nrows)) (1/ncols) (1/nrows)]);
%axes('position',[(x-1)*(1/ncols) 1-(y*(1/nrows)) (1/ncols) (1/nrows)]);
% TO DO:  fix this so it works in dual-representation mode

      lego(z{fi},P);
      if ~isempty(titles)
        zl=zlim;
        text(1,1,zl(2)*1.6,titles{fi},'fontsize',P.fontsize,'interpreter','none','verticalalignment','top');
      end
    end
    fi=fi+1;
end,end,set(gcf,'color',[1 1 1]);

