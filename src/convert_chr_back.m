function L = convert_chr_back(n,P)
% convert_chr_back(chromosome_list)
%
% converts numbers 1-24 to chromosome identifiers (1-->chr1, 23-->chrX, 24-->chrY)
%
% Mike Lawrence 2009-05-19

if ~exist('P','var'), P=[]; end
if ischar(P)
  tmp = P;
  P = [];
  P.build = tmp;
end
P = impose_default_value(P,'build','');

% CONVERSION METHODS AVAILABLE:
% if "build" is not specified, asssume human.
% if "build" is specified, but ReferenceInfoObj has not been initialized, use heuristics defined here.
% if "build" is specified AND ReferenceInfoObj has been initialized, use it to do the conversion.

RIO_available = false;

if isempty(P.build)
  assumed_human = true;
  ct = 24;
else
  try
    ct = ReferenceInfoObj.getMaxNum(P.build);
    RIO_available = true;
  catch me
    assumed_human = false;
    ct = get_chrcount(P.build);
  end
end

% CONVERSION

if isempty(n)
  L = [];
  return
end

if iscell(n) && ~strncmp(n{1},'chr',3)
  n = str2double(n);
end

if ~isnumeric(n)
  L = n;
  return
end

L = cell(size(n));

idx = find(n<0 | n>ct);
if ~isempty(idx)
  fprintf('WARNING: chromosome number(s) out of range for this genome:\n');
  disp(unique(n(idx)));
  keyboard
end

if ~RIO_available
  % legacy behavior
  
  for i=1:length(n)
    L{i} = ['chr' num2str(n(i))];
  end
  
  L = regexprep(L,'chr0','chrM');
  idx = grep(['^chr' num2str(ct-1) '|chr' num2str(ct) '$'],L,1);
  if ~isempty(idx)
    if assumed_human
      fprintf('convert_chr: assuming human for chrX/chrY\n');
    end
    L = regexprep(L,['chr' num2str(ct-1)],'chrX');
    L = regexprep(L,['chr' num2str(ct)],'chrY');
  end

else
  % use ReferenceInfoObj
  
  for i=1:length(n)
    zz = num2str(n(i));
    try
      if strcmp(ReferenceInfoObj.getUse(n(i),P.build),'D')
        zz = decell(ReferenceInfoObj.getChrom(n(i),P.build));
      end
    catch me
    end
    if ~strncmp(zz,'chr',3), zz = ['chr' zz]; end
    L{i} = zz;
  end
end



