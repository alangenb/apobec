function lego(z,P)
% lego(vals96,names96)   names in format "A in C_T ->G"
%
% legacy input:
%        8x12 matrix of heights for the "LEGO" plot
%        tries to convert other formats

if ~exist('P','var'), P=[]; end

if iscellstr(P)   % can specify catnames as second argument for factor2lego
  tmp=P;
  P=[];
  P.catnames = tmp;
end

if iscell(z), legos(z,P); return; end
if isstruct(z) && isfield(z,'mut'), z=z.mut; end
if ischar(z) || (isstruct(z) && (isfield(z,'context65')||(isfield(z,'chr')&(isfield(z,'start')|isfield(z,'st')|isfield(z,'pos'))))), legomaf(z,P); return; end
if length(size(z))==2 & all(size(z)==[1 96]), z = factor2lego(z,P); end
if length(size(z))==2 & all(size(z)==[96 1]), z = factor2lego(z',P); end
if length(size(z))==4 & all(size(z)==[4 4 4 4] | size(z)==[4 2 4 4]), z = matrix2lego(z); end
if ~(length(size(z))==2 & all(size(z)==[8 12])), error('unknown format'); end

if ischar(P) && ismember(P,{'counts','rates','exome','genome'})
  if ~strcmp(P,'counts'), error('only "counts" supported'); end
  P=[];
end
 
P = impose_default_value(P,'log',false);
P = impose_default_value(P,'normalize',true);
P = impose_default_value(P,'imagesc',false);  % if true, displays as imagesc instead of bar3_with_colors
P = impose_default_value(P,'zmax',[]);
if isfield(P,'show_sanger_plot')&P.show_sanger_plot, P.add_sanger_plot=true;end
if isfield(P,'dualrep')&P.dualrep, P.add_sanger_plot=true;end
P = impose_default_value(P,'add_sanger_plot',false);
P = impose_default_value(P,'normalize',true);

if P.log, z=log(z); end

if P.normalize, z = z/sum(z(:)); end

if P.imagesc
  if ~isempty(P.zmax), z = min(z,P.zmax); end
  imagesc(z);
  pp = {'linewidth',2,'color',[0 0 0]};
%
% blank lines for code alignment with lego5.m
%
  line([4.5 4.5],ylim,pp{:});line([8.5 8.5],ylim,pp{:});line(xlim,[4.5 4.5],pp{:});
else
  c = get_LEGO_colors(P);
%  z3 = repmat(z,[1 1 3]); c(isnan(z3))=1; % white   %%% for some reason this is causing bizarre color scrambling in context of legos.m
  bar3_with_colors(z,c);
  if ~isempty(P.zmax)
    zlim([0 P.zmax]);
  else
    zlim([0 1.1*max(z(:))]);
  end
  % NOTE: sometimes in Mac Xwin "Xquartz", the plot disappears if zmax drops below ~0.4e-6
  zl = zlim;
  if zl(2)<0.4e-6
    fprintf('Warning: plot sometimes disappears when zmax<0.4e-6.  If this happens, try typing zlim([0 0.4e-6]).\n');
  end
end
set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);

if P.add_sanger_plot
  set(gca,'position',[0.13 -0.03 0.775 0.815]);  % move main plot slightly down
  h = axes('position',[0.05 0.77 0.88 0.22]);
  set(h,'visible','off','xtick',[],'ytick',[]);
  sz = lego_to_sanger(z);
  sc = get_sanger_colors();
  sz = sz/max(sz);
  for i=1:96
    h = 0.9*sz(i); if h==0, continue; end
    y = 0.02;
    x = 0.022+(0.94/96)*i;
    w = (0.94/96)*0.8;
    rectangle('position',[x y w h],'facecolor',sc(i,:),'linestyle','none');
  end
  rectangle('position',[0.02 0.02 0.96 0.96],'linewidth',1,'edgecolor',[0.8 0.8 0.8]); % bounding box
end




