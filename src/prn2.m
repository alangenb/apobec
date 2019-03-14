function output = prn(varargin)
% prn(...)
% 
% same as pr('#',...)
% but easier to type

if nargout==0
  pr('#',varargin{:});
elseif nargout==1
  output = pr('#',varargin{:});
else
  error('too many output variables');
end


