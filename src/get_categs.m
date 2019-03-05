function [Z,C] = get_categs(categdir,build)
% [Z,C] = get_categs(categdir)
%
% Z = list of category names
% C = cell(24,1) of category data arrays [only loaded if requested]

if ~exist('build','var'), build='hg19'; end

if categdir(1)=='/'
  dr = categdir;
else
  dr = ['/cga/tcga-gsc/home/lawrence/db/' build '/' categdir];
end

Z = load_struct([dr '/categs.txt']);

if nargout>=2
  fprintf('Loading categ data from %s\n',dr);
  C = cell(24,1);
  for i=1:24, fprintf('chr%d ',i);
    tmp = load([dr '/chr' num2str(i) '.mat']);
    f = fieldnames(tmp);
    C{i} = getfield(tmp,f{1});
  end,fprintf('\n');
end


