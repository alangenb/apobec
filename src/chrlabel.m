function varargout = chrlabel(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = convert_chr_back(varargin{:});
else
  [varargout{1}] = convert_chr_back(varargin{:});
end
