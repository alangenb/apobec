function X = make_input_for_nmf(m,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'coding_only',true);

if ischar(m)
  if exist(m,'dir')
    m = load2(m);
  elseif exist(m,'file')
    m = load_struct(m);
  else
    error('Not found: %s',m);
  end
end

% keep only coding
if P.coding_only
  fprintf('Subsetting to coding mutations\n');
  if isfield(m,'type')
    m = reorder_struct(m,grepmi('missense|nonsense|splice|synon|silent|ins|del',m.type));
  elseif isfield(m,'is_coding')
    m = reorder_struct(m,m.is_coding);
  else
    error('Need m.type to subset to coding only');
  end
end

if ~isfield(m,'context65')
  m.context65 = get_context65(m.chr,m.pos);
end
m = make_numeric(m,'context65');
orig_list = generate_categ_context65_names();

if ~isfield(m,'newbase_idx')
  if isfield(m,'alt_idx')
    m = rename_field(m,'alt_idx','newbase_idx');
  elseif isfield(m,'alt')
    if isnumeric(m.alt)
      m.newbase_idx = m.alt;
    elseif iscellstr(m.alt)
      m.newbase_idx = listmap(upper(m.alt),{'A','C','G','T'});
    end
  else
    m.newbase_idx = listmap(upper(m.newbase),{'A','C','G','T'});
  end
end

if isfield(m,'classification')
  m.is_snp = strcmp('SNP',m.classification);
else
  m.is_snp = ~isnan(m.newbase_idx);
end

% patients
[tmp pati m.pat_idx] = unique(m.patient);
pat = reorder_struct(keep_fields_if_exist(m,{'patient','ttype','dataset','source','repo','ttype0'}),pati);
pat = rename_field(pat,'patient','name'); pat=order_field_first(pat,'name');

% SNPs
midx = (m.is_snp & m.pat_idx>=1 & m.pat_idx<=slength(pat) & m.context65>=1 & m.context65<=65 & m.newbase_idx>=1 & m.newbase_idx<=4);
pat.nsnp = hist3d_fast(m.pat_idx(midx),m.context65(midx),m.newbase_idx(midx),1,slength(pat),1,65,1,4);
pat.nchan = nan(slength(pat),1);

% convert 65x4 to 192
pat.nchan = zeros(slength(pat),192);
chan=[]; chan.name = cell(192,1);
base='ACGT';
i=1;
for from=1:4
  for to=1:4
    for left=1:4
      for right=1:4
        if from==to, continue; end
        orig_name = [base(from) ' in ' base(left) '_' base(right)];        
        chan.name{i,1} = [base(left) ' (' base(from) '->' base(to) ') ' base(right)];
        cidx = find(strcmp(orig_name,orig_list.name)); if isempty(cidx), error('what?'); end
        pat.nchan(:,i) = pat.nsnp(:,cidx,to);
        i=i+1;
end,end,end,end

% indels+DNPs
if isfield(m,'classification')
  midx = find(strcmp('DEL',m.classification));
  pat.ndel = histc(m.pat_idx(midx),1:slength(pat));
  midx = find(strcmp('INS',m.classification));
  pat.nins = histc(m.pat_idx(midx),1:slength(pat));
  midx = find(strcmp('DNP',m.classification));
  pat.ndnp = histc(m.pat_idx(midx),1:slength(pat));
end

X=[];
X.pat = rmfield(pat,'nsnp');
X.chan = chan;

X = collapse_nmf_input_192_to_96(X);
X.chan.catnames = regexprep(X.chan.name,'^(.) \((.)->(.)\) (.)$','$2 in $1_$4 ->$3');

% check for patients with zero mutations
totmut = sum(X.pat.nchan,2);
idx = find(totmut==0);
if ~isempty(idx)
  fprintf('WARNING: removing %d/%d patients that have zero mutations.\n',length(idx),slength(X.pat));
  X.pat = reorder_struct_exclude(X.pat,idx);
end

