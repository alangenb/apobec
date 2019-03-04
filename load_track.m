function [Z,C] = load_track(varargin)

if nargout<2
 Z = get_categs(varargin{:});
elseif nargout==2
 [Z,C] = get_categs(varargin{:});
else
  error('too many outputs requested');
end

