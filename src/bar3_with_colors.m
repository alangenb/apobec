function bar3_with_colors(varargin)
% bar3_with_colors(varargin)
%
% all parameters but the last one are passed unchanged to bar3()
% last parameter should be matrix C, the same x*y size as the matrix being plotted, with size(C,3)==3,
% with page1,2,3=RGB for each matrix element
%
% future modification plans:
%    - individual control of the individual faces of each bar (encoded in the 6*nrows)
%
% Mike Lawrence 2010-11-08

h = bar3(varargin{1:end-1});
ncols = length(h);
nrows = size(get(h(1),'zdata'),1)/6;

C = varargin{end};
if size(C,1)~=nrows || size(C,2)~=ncols || size(C,3)~=3, error('C is wrong dimensions'); end

allcolors = zeros(ncols*nrows,3);
idx=1;for y=1:nrows, for x=1:ncols, allcolors(idx,:)=C(y,x,:); idx=idx+1;end;end
[u ui uj] = unique(allcolors,'rows');

colormap(gca, u);
for x=1:ncols
  xdata=get(h(x),'xdata');
  ydata=get(h(x),'ydata');
  zdata=get(h(x),'zdata');
  cdata=get(h(x),'cdata');
  for y=1:nrows
    cdata(y*6-5:y*6,:) = uj((y-1)*ncols+x);
    if isnan(varargin{1}(y,x))  % make height=NaN bars invisible
      xdata(y*6-5:y*6,:) = nan;
      ydata(y*6-5:y*6,:) = nan;
      zdata(y*6-5:y*6,:) = nan;
    end
  end
  set(h(x),'xdata',xdata,'ydata',ydata,'zdata',zdata,'cdata',cdata);
end

