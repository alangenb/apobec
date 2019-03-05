function varargout = mmw(varargin)
if nargout==0
  map_mutations_to_targets_fast_wrapper(varargin{:});
elseif nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = map_mutations_to_targets_fast_wrapper(varargin{:});
else
  [varargout{1}] = map_mutations_to_targets_fast_wrapper(varargin{:});
end
