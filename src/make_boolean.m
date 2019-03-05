function varargout = make_boolean(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = make_logical(varargin{:});
else
  [varargout{1}] = make_logical(varargin{:});
end
