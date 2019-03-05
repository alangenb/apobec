function varargout = write_textfile(varargin)
if nargout==0
  save_textfile(varargin{:});
elseif nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = save_textfile(varargin{:});
else
  [varargout{1}] = save_textfile(varargin{:});
end
