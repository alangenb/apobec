function varargout = load_categs(varargin)
if nargout==0
  get_categs(varargin{:});
elseif nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = get_categs(varargin{:});
else
  [varargout{1}] = get_categs(varargin{:});
end
