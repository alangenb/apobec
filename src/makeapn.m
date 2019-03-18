function varargout = makeapn(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = make_all_possible_numeric(varargin{:});
else
  [varargout{1}] = make_all_possible_numeric(varargin{:});
end
