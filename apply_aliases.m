function list = apply_aliases(list, aliases, varargin)
%
% apply_aliases(list, aliases, parameters)
%
% Renames any names in list that match aliases.old to 
%   the corresponding names in aliases.new
%
% If 'ignore_extensions' appears among parameters,
%    then spliceform identifiers (.*) are ignored and preserved,
% 
% If 'case_insensitive' appears among parameters,
%    then case is ignored in making the comparisons
%
% Mike Lawrence 2008-05-22
%

ignore_extensions = false;
case_insensitive = false;

for i=1:length(varargin)
  if strcmpi(varargin{i}, 'ignore_extensions')
    ignore_extensions = true;
  elseif strcmpi(varargin{i}, 'case_insensitive')
    case_insensitive = true;
  else
    error('Unrecognized option %s', varargin{i});
  end
end

require_fields(aliases, {'old';'new'});

if length(aliases.old) ~= length(aliases.new)
  error('old and new namelists must be of equal length.\n');
end

if ignore_extensions
  extensions = regexp(list,'(\..*)$','match');
  list = regexprep(list,'(\..*)','');
end

for i=1:length(aliases.old)
  if case_insensitive
    match = strcmpi(list, aliases.old{i});
  else
    match = strcmp(list, aliases.old{i});
  end
  list(match) = repmat(aliases.new(i), sum(match), 1);
end

if ignore_extensions    % put them back
  for i=1:length(list)
    for j=1:length(extensions{i})
      list{i} = [list{i} extensions{i}{j}];
    end
  end
end
