function print_to_file(filename,varargin)

if nargin == 1
  [dev res] = interpret_print_filename(filename);
else 
  dev = interpret_print_filename(filename);
  res = varargin{1};
end 
%keyboard

if nargin == 1
  [dev res] = interpret_print_filename(filename);
else
  dev = interpret_print_filename(filename);
  res = varargin{1};
end

fprintf('Outputting figure to %s\n', filename);
print(['-d' dev], ['-r' num2str(res)], filename);

end
