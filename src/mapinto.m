function X = mapinto(X,Y,xkey,ykey,yfields,yfields_rename)
% X = mapinto(X,Y,xkey[,ykey[,yfields[,yfields_rename]]])
% Mike Lawrence 2011-2012

if nargin==3
  % special-case usage:
  % X = mapinto(X,Y,key)
elseif nargin==4
  if ischar(xkey) && ischar(ykey) && isfield(X,xkey) && isfield(Y,ykey)
    % special-case usage:
    % X = mapinto(X,Y,xkey,ykey)
  elseif ischar(xkey) && ischar(ykey) && isfield(X,xkey) && ~isfield(Y,ykey)
    % special-case usage:
    % X = mapinto(X,Y,key,yfield)
    yfields = ykey;
    clear ykey;
  elseif iscell(ykey)
    % special-case usage:
    % X = mapinto(X,Y,key,yfields)
    yfields = ykey;
    clear ykey;
  end
elseif nargin==5
  if ischar(xkey) && ischar(ykey) && isfield(X,xkey) && isfield(Y,ykey)
    % special-case usage:
    % X = mapinto(X,Y,xkey,ykey,yfields)
  elseif iscell(ykey) || (ischar(yfields) && ~isfield(Y,yfields))
    % special-case usage:
    % X = mapinto(X,Y,key,yfields,yfields_rename)
    yfields_rename = yfields;
    yfields = ykey;
    clear ykey;
  end
elseif nargin==6
  % full-usage case: X = mapinto(X,Y,xkey,ykey,yfields,yfields_rename)
else
  error('usage: X = mapinto(X,Y,xkey[,ykey[,yfields[,yfields_rename]]])');
end

if ~exist('ykey','var')
  ykey = xkey;
end

if ~exist('yfields','var')
  yfields = fieldnames(Y);
end

if exist('yfields_rename','var')
  if length(yfields)~=length(yfields_rename), error('length(yfields)~=length(yfields_rename)'); end
else
  yfields_rename = yfields;
end

if ~isstruct(X) || ~isstruct(Y)
  error('X and Y should be structs');
end

if ~ischar(xkey) || ~ischar(ykey)
  error('xkey and ykey should be strings');
end

if ischar(yfields)
  yfields = {yfields};
end
if ~iscell(yfields)
  error('yfields should be a string or cell array of strings');
end

xx = getfield(X,xkey);
yy = getfield(Y,ykey);
idx = listmap(xx,yy);
ii = find(~isnan(idx));
if isempty(ii), fprintf('Warning: mapinto matched no entries\n'); end

for i=1:length(yfields), fld=yfields{i};,fld2=yfields_rename{i};
  if strcmp(fld,ykey), continue; end  % don't overwrite the key
  yy = getfield(Y,fld);
  if isfield(X,fld2)
    xx = getfield(X,fld2);
    if iscell(xx) && isnumeric(yy)
      fprintf('Need to resolve types for %s: converting numeric to cellstr\n',fld);
      yy=num2cellstr(yy);
    elseif isnumeric(xx) && iscell(yy)
      fprintf('Need to resolve types for %s: converting numeric to cellstr\n',fld);
      xx=num2cellstr(xx);
    end
    xs = size(xx); ys = size(yy);
    if length(xs)~=length(ys) || any(xs(2:end)~=ys(2:end))
      whos xx yy
      error(['Incompatible inputs in field "' fld '"']);
    end
    xx(ii,:,:,:,:,:,:,:,:,:) = yy(idx(ii),:,:,:,:,:,:,:,:,:);  % don't put overwrite previous values with NaN where missing
  else
    xx = nansub(yy,idx);
  end
  X = setfield(X,fld2,xx);
end

