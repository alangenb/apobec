function R = make_all_possible_numeric(R,threshold)
% R = make_all_possible_numeric(R,threshold)
%
% For each field, tries to make it numeric.
% if at least <threshold> records become non-NaN, then switches to the numeric version.
%
% Mike Lawrence 2010-04-20

if ~exist('threshold','var'), threshold = 0.99; end
if threshold>1
  fprintf('Converting threshold=%d to threshold=%d%%\n',threshold,threshold/100);
  threshold = threshold / 100;
end

len = slength(R);

f = fieldnames(R);
for i=1:length(f)
  fprintf('%-20s -->  ', f{i});
  x = getfield(R,f{i});
  if isnumeric(x) || islogical(x)
    fprintf('Already numeric');
  else try
    wasnan = strcmpi(x,'NaN') | strcmpi(x,'Inf') | cellfun('isempty',x);
    xn = str2doubleq_wrapper(x);
    tonum = (wasnan|~isnan(xn));
    frac = sum(tonum)/len;
    if frac==0, fprintf('All non-numeric  ');
    elseif frac==1, fprintf('All numeric  ');
    else fprintf('%d/%d numeric (%0.1f%%) ',sum(tonum),len,100*frac);
    end
    if frac>=threshold
      fprintf('--> converted');
      R = setfield(R,f{i},xn);
    else
      fprintf('--> left as is');
    end
  catch me
    fprintf('Error: %s',me.message);
  end, end
  fprintf('\n');
end
