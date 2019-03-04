function x = get_LEGO_colors(P)

if ~exist('P','var'),P=[];end

P = impose_default_value(P,'use_lighter_blue',false);

if P.use_lighter_blue % for CMYK-space print figures
  colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0.2 0.5 0.95;0.7 0.3 0.75];
else
  colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7];
end

t = [repmat([3 3 3 3 2 2 2 2 1 1 1 1],4,1);repmat([6 6 6 6 5 5 5 5 4 4 4 4],4,1)];
r = colors(:,1); g = colors(:,2); b = colors(:,3);
x = cat(3,r(t),g(t),b(t));



