function [z c] = factor2lego(h,P)
% [z c] = factor2lego(h,catnames)
% [z c] = factor2lego(h,P) with P.catnames
%
% Helper function for LEGO-plot visualization of NMF factors
%
% GIVEN "h", a set of factor definitions from NMF as performed in mutation_spectra_pca_plot.m
%          each row = a factor
%          columns = the 96 categories of mutations (after strand collapse)
% and "catnames"   e.g. "A in C_T ->G"
%         = tells what order the 96 categories are in.  (default order used if not provided)
%
% REORDERS "h" to yield:
%       z = a set of 8x12 heights for LEGO plot
%
% RETURNS also:
%       c = a set of 8x12 colors for LEGO plot
%
% Note: if "h" has multiple rows, then "z" will have multiple pages
%       (but c will not be provided in duplicate)
%
% Mike Lawrence 2012-05-22

if ~exist('P','var')
  P=[];
elseif iscellstr(P)
  tmp=P; 
  P=[];
  P.catnames = tmp;
elseif isstruct(P) || isempty(P)
  % OK
elseif ischar(P) && ismember(P,{'counts','rates','exome','genome'});
  if ~strcmp(P,'counts'), error('only "counts" supported'); end
  P=[];
else
  error('invalid second parameter');
end

P = impose_default_value(P,'catnames',generate_lego_96names());

if size(h,2)~=96, error('h should have 96 columns'); end

classes = {'C->G','C->A','C->T','A->T','A->C','A->G'};
bases = 'TCAG';

z = zeros(8,12,size(h,1));
cx=1; cy=1;
for ci=1:6
  for left=1:4
    for right=1:4
      name = [classes{ci}(1) ' in ' bases(left) '_' bases(right) ' ' classes{ci}(2:end)];
      cidx = find(strcmp(P.catnames,name));
      if length(cidx)~=1
        name = [bases(left) ' (' classes{ci} ') ' bases(right)];
        cidx = find(strcmp(P.catnames,name));
        if length(cidx)~=1
          error('Problem with P.catnames');
        end
      end
      y = cy+left-1;
      x = cx+right-1;
      z(y,x,:) = h(:,cidx);
  end,end
  cx=cx+4;
  if cx>9, cx=1;cy=cy+4; end
end

if nargout>=2
  c = get_LEGO_colors();
end

















