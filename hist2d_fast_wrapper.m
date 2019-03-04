function h = hist2d_fast_wrapper(x,y,xmin,xmax,ymin,ymax)

if nargin~=2 && nargin~=6
  error('needs 2 or 6 inputs');
end

x = real(double(x));
y = real(double(y));

if any(isnan(x)) || any(isnan(y))
  if nargin==2
    error('Can''t handle NaNs properly without specifying xmin,xmax,ymin,ymax');
  elseif nargin==6
    xmin = real(double(xmin));
    xmax = real(double(xmax));
    ymin = real(double(ymin));
    ymax = real(double(ymax));
    x(isnan(x)) = xmin-1;
    y(isnan(y)) = ymin-1;
  end
end

if nargin==6
  h = hist2d_fast(x,y,xmin,xmax,ymin,ymax);
elseif nargin==2
  h = hist2d_fast(x,y);
end
