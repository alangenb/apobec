function varargout = as_col(varargin)
if nargout==0
  as_column(varargin{:});
elseif nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = as_column(varargin{:});
else
  [varargout{1}] = as_column(varargin{:});
end
