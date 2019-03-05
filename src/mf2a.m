function varargout = mf2a(varargin)
if nargout==0
  move_field_to_after(varargin{:});
elseif nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = move_field_to_after(varargin{:});
else
  [varargout{1}] = move_field_to_after(varargin{:});
end
