function a = get_context(chr,pos,contextdir,P)

if nargin<3, error('need at least 3 inputs'); end

if isempty(chr)||isempty(pos), fprintf('empty query\n'); a=[]; return; end
 
if ischar(chr)
  fprintf('Assuming %s is contextdir\n',chr);
  tmp=chr;
  chr=pos;
  pos=contextdir;
  contextdir=tmp;
end

if exist('P','var') && ischar(P)
  fprintf('get_context fourth argument was char:  assuming this means fileext');
  tmp = P;
  P=[];
  P.fileext = tmp;
end

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'build','');
P = impose_default_value(P,'fileext','.mat');

if ~isnumeric(chr), chr = convert_chr(chr,P.build); end
if ~isnumeric(pos), pos = str2double(pos); end

if length(chr)==1 && length(pos)>1
  chr = repmat(chr,length(pos),1);
end

if length(chr) ~= length(pos), error('chr and pos should be vectors of the same length'); end

% CHOOSE WHAT FORMAT TO READ

% Two possible formats to pull from:  FWB and MAT.
% For moderate-size query sets, FWB is faster, but it might not be available
% For very large query sets, MAT is faster, but it might not be available
% Find out which are available and choose accordingly.
query_size_threshold = 1e4;

fwb = [contextdir '/all.fwb'];
fwb_available = exist(fwb,'file');

mat_available = true;
for c=1:24
  mat = [contextdir '/chr' num2str(c) P.fileext];
  if ~exist(mat,'file')
    mat_available = false;
    break;
  end
end

format_to_use = '';
if ~fwb_available && ~mat_available
  error('Can''t find either FWB or MAT file in the context directory!');
elseif fwb_available && ~mat_available
  format_to_use = 'fwb';
elseif mat_available && ~fwb_available
  format_to_use = 'mat';
else % both available
  if length(chr)<query_size_threshold
    format_to_use = 'fwb';
  else
    format_to_use = 'mat';
  end
end

% RETRIEVE DATA

a = nan(length(chr),1);

fprintf('Getting context: ');

if strcmp(format_to_use,'fwb')
  fwb = [contextdir '/all.fwb'];
  try
    if isempty(P.build)
      a = get_from_fwb(fwb,chr,pos);
    else
      a = get_from_fwb(fwb,chr,pos,P.build);
    end
  catch me
    if ~mat_available
      disp(me);
      disp(me.message);
      error('Error in get_from_fwb');
    else
      format_to_use = 'mat';
    end
  end
end


if strcmp(format_to_use,'mat')
  for c=1:24
    idx = find(chr==c);
    if isempty(idx), continue; end
    fprintf('chr%d ',c);
    x = load([contextdir '/chr' num2str(c) P.fileext]);
    f = fieldnames(x);
    x = getfield(x,f{1});
    maxlen = length(x);
    idx1 = idx(pos(idx)<=maxlen);
    idx2 = idx(pos(idx)>maxlen);
    a(idx1) = x(pos(idx1));
    a(idx2) = nan;
  end
end

fprintf('\n');


