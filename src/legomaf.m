function [z c] = legomaf(X,P)
% legomaf(X,P)

trimer_coverage_file = '/cga/tcga-gsc/home/lawrence/db/hg19/context65/breakdowns.mat';

if exist('P','var') && ischar(P)
  tmp=P;
  P=[];
  P.coverage = tmp;
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'display','counts');
P = impose_default_value(P,'coverage','exome');

%if grepmi('exome',{P.coverage}), N = average_exome_coverage();
%elseif grepmi('genome|wgs',{P.coverage}), N = average_genome_coverage();
load(trimer_coverage_file,'K');

if grepmi('^(flat|counts|n)$',{P.coverage})||strcmp('counts',P.display), P.display='counts'; N = 1e6*ones(64,1); fprintf('Showing raw counts (no normalization to territory or coverage)\n');
elseif grepmi('exome|rates',{P.coverage}), load(trimer_coverage_file,'K'); N = K.Nx(1:64); fprintf('Normalizing by exome territory\n');
elseif grepmi('genome|wgs',{P.coverage}), load(trimer_coverage_file,'K'); N = K.Ng(1:64); fprintf('Normalizing by genome territory\n');
else
  error('don''t know what coverage scheme to use');
end

if ischar(X), X = load_struct(X); end

if isfield(X,'classification'), X = reorder_struct(X,strcmp('SNP',X.classification)); end
if isfield(X,'Variant_Type'), X = reorder_struct(X,strcmp('SNP',X.Variant_Type)); end
if isfield(X,'is_indel'), X = make_boolean(X, 'is_indel'); X = reorder_struct(X,~X.is_indel); end

%%%%%%%%%

if ~isfield(X,'newbase_idx') && ~isfield(X,'newbase') && isfield(X,'alt')
  if isnumeric(X.alt), X.newbase_idx = X.alt;
  elseif iscellstr(X.alt), X.newbase = X.alt;
  end
end

if ~isfield(X,'newbase_idx')
  if ~isfield(X,'newbase')
    X.newbase = find_newbase(X);
  end
  X.newbase_idx = listmap(X.newbase,{'A','C','G','T'});
end

if ~isfield(X,'context65')
  if ~isfield(X,'pos') && isfield(X,'start'), X.pos=X.start; end
  if ~isfield(X,'pos') && isfield(X,'st'), X.pos=X.st; end
  X.context65 = get_context65(X.chr,X.pos);
end

X = make_numeric(X,{'context65','newbase_idx'});





midx = find(X.context65>=1 & X.context65<=64 & X.newbase_idx>=1 & X.newbase_idx<=4);
try
  n = hist2d_fast_wrapper(X.context65(midx),X.newbase_idx(midx),1,64,1,4);
catch
  n = hist2d(X.context65(midx),X.newbase_idx(midx),1:64,1:4);
end

%%%%%%%%

if isfield(X,'patient'), npat = length(unique(X.patient));
elseif isfield(X,'pat_idx'), npat = length(unique(X.pat_idx));
elseif isfield(X,'patient_idx'), npat = length(unique(X.patient_idx));
elseif isfield(X,'Tumor_Sample_Barcode'), npat = length(unique(X.Tumor_Sample_Barcode));
else npat=1; end

N=N*npat;

Nn = [N n];

Nn = collapse_Nn_64_by_strand(Nn);

% convert to 96-row format:
% 17-32 >A
% 1-16 >C
% 1-16 >G
% 17-32 >G
% 1-16 >T
% 17-32 >T
N = [Nn(17:32,1);Nn(1:16,1);Nn(1:16,1);Nn(17:32,1);Nn(1:16,1);Nn(17:32,1)];
n = [Nn(17:32,2);Nn(1:16,3);Nn(1:16,4);Nn(17:32,4);Nn(1:16,5);Nn(17:32,5)];

if strcmp(P.display,'counts'), data = n;
elseif strcmp(P.display,'rates'), data = n./N;
else error('Unknown P.display %s',P.display); end

% Lego plot
z = factor2lego(data');
lego(z,P);

if nargout==0
  clear z c
end
