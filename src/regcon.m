function out = regcon(input,regex)
%find instances of elements in cellstr containg regex
%usage: output = regcon(input,regex);
if ~iscellstr(input) error('Error: input must be cellst'); end
if ~ischar(regex) error('Error: regex must be a string'); end
out = cell2mat(cellfun(@(x) length(x)~=0, regexp(input,regex,'match','once'), 'uni', 0));
end
