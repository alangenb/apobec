function G = get_refseq_introns_and_exons(build,splicesiteflank,P)
% G = get_refseq_introns_and_exons
%
% returns list of RefSeq genes and their exons and introns
% notes:
%   (1) exons here means "coding sequences"; ignores UTRs
%   (2) introns exclude regions that are exonic in *any* gene
%   (3) G.name is not necessarily unique, because some genes appear on more than one chromosome
%
% Mike Lawrence 2010-06-19

if ~exist('build','var'), build = 'hg19gencode'; fprintf('Assuming build %s\n',build); end

if exist('splicesiteflank','var') && isstruct(splicesiteflank) && ~exist('P','var')
  P = splicesiteflank;
  clear splicesiteflank;
end

if ~exist('splicesiteflank','var'), splicesiteflank = 0; end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'genes_to_ignore',[]);
P = impose_default_value(P,'add_promoters',false);
P = impose_default_value(P,'imputed_promoter_length',3000);
P = impose_default_value(P,'skip_longname_annotation',false);
P = impose_default_value(P,'load_conservation29',false);
P = impose_default_value(P,'load_conservation3',false);

premadefile=['/cga/tcga-gsc/home/lawrence/db/load_genes/' build '_splicesiteflank' num2str(splicesiteflank) '.mat'];
if exist(premadefile,'file')
  load(premadefile,'G');

else

  R = load_refseq(build);
  if strcmp(build,'hg19gencode'), build='hg19'; end
  if grepmi('gencode|gc19', {build}), build = 'hg19'; end
  
  if ~isempty(P.genes_to_ignore)
    ridx = find(ismember(R.gene,P.genes_to_ignore));
    if ~isempty(ridx)
      fprintf('Omitting the following genes:\n'); pr(unique(R.gene(ridx)));
      R = reorder_struct_exclude(R,ridx);
    end
  end
  
  R = reorder_struct(R, ~strcmp(R.chr, 'chrM'));
  R.chr = convert_chr(R.chr,build);
  R = reorder_struct(R,~isnan(R.chr));
  R = reorder_struct(R,R.code_len>0);
  
  nchr = get_chrcount(build);
  
  fprintf('Collapsing transcripts for each gene... ');
  G = cell(nchr,1);
  for chr=1:nchr, fprintf('chr%d ',chr);
    Rc = reorder_struct(R,R.chr==chr);
    [ug ugi ugj] = unique(Rc.gene);
    % (1) for each gene, get list of coding exons (uniqued and collapsed)
    G{chr}.name = ug; G{chr}.chr = chr*ones(length(ug),1);
    G{chr}.strand = cell(length(ug),1);
    G{chr}.tx_start = nan(length(ug),1); G{chr}.tx_end = nan(length(ug),1);
    G{chr}.code_start = nan(length(ug),1); G{chr}.code_end = nan(length(ug),1);
    G{chr}.n_exons = nan(length(ug),1); G{chr}.exon_starts = cell(length(ug),1); G{chr}.exon_ends = cell(length(ug),1);
    G{chr}.n_introns = nan(length(ug),1); G{chr}.intron_starts = cell(length(ug),1); G{chr}.intron_ends = cell(length(ug),1);  
    for gi=1:length(ug), idx = find(ugj==gi);
      G{chr}.strand{gi} = concat(unique(Rc.strand(idx)),'/');
      G{chr}.tx_start(gi) = min(Rc.tx_start(idx)); G{chr}.tx_end(gi) = max(Rc.tx_end(idx));
      E = cell(length(idx),1);
      for i=1:length(idx), j=idx(i);
        E{i} = [max(Rc.code_start(j), Rc.exon_starts{j}) min(Rc.code_end(j), Rc.exon_ends{j})];
        E{i}(E{i}(:,2)<E{i}(:,1),:) = [];  % remove noncoding exons
        E{i}(:,1) = E{i}(:,1) - splicesiteflank;
        E{i}(:,2) = E{i}(:,2) + splicesiteflank;
      end
      E = unique(cat(1,E{:}),'rows'); % (sorts by start,end)
      while(true)   % collapse overlapping/adjacent exons
        idx = find(E(2:end,1)<=E(1:end-1,2)+1,1); if isempty(idx), break; end
        E(idx,2) = max(E(idx:idx+1,2)); E(idx+1,:) = [];
      end
      G{chr}.code_start(gi) = E(1,1); G{chr}.code_end(gi) = E(end,2);
      G{chr}.n_exons(gi) = size(E,1); G{chr}.exon_starts{gi} = E(:,1); G{chr}.exon_ends{gi} = E(:,2);
    end
    % (2) for each gene,
    %    get list of introns as follows:
    %    -- take interval from coding start to coding end
    %    -- remove intervals that are exons of this or any other gene
    for gi=1:length(ug)
      idx = find(G{chr}.code_start<=G{chr}.code_end(gi) & G{chr}.code_end>=G{chr}.code_start(gi));
      I = [G{chr}.code_start(gi) G{chr}.code_end(gi)];
      for i=1:length(idx), j=idx(i);
        for e=1:length(G{chr}.exon_starts{j})
          est = G{chr}.exon_starts{j}(e); een = G{chr}.exon_ends{j}(e);
          ii = find(I(:,1)<=een & I(:,2)>=est);
          for k=1:size(ii,1), iik = ii(k);
            if I(iik,1)>=est && I(iik,2)<=een    % whole intron is exonic
              I(iik,:) = [nan nan];
            elseif I(iik,1)<est && I(iik,2)>een   % intron is bisected by an exon
              I(end+1,:) = [een+1 I(iik,2)];
              I(iik,2) = est-1;
            elseif I(iik,2)>een     % intron needs to be trimmed after exon
              I(iik,1) = een+1;
            else                   % intron needs to be trimmed before exon
              I(iik,2) = est-1;
            end,end,end,end
            I = sortrows(I(~isnan(I(:,1)),:)); G{chr}.n_introns(gi) = size(I,1);
            G{chr}.intron_starts{gi} = I(:,1); G{chr}.intron_ends{gi} = I(:,2);
    end
    
  end, fprintf('\n');  % next chromosome
  
  G = concat_structs(G);
  if ~P.skip_longname_annotation,
    G.longname = get_longnames(G.name);
  else
    G.longname = repmat({''}, slength(G), 1);
  end
  G = order_fields_first(G,{'name','longname'});
  
  for i=1:slength(G), G.tot_intron_len(i,1) = sum(G.intron_ends{i}-G.intron_starts{i}+1); end
  G.footprint = G.code_end-G.code_start+1;
  G.footprint_including_UTRs = G.tx_end-G.tx_start+1;
  G.tot_exon_len = G.footprint-G.tot_intron_len;
  G.coding_len = nan(slength(G),1);
  for i=1:slength(G), G.coding_len(i) = sum(G.exon_ends{i}-G.exon_starts{i}+1); end

  save(premadefile,'G');

end
  
% add promoters if requested

G.gene_start = G.tx_start;
G.gene_end = G.tx_end;
if P.add_promoters
  idx = grep('+',G.strand,1);
  G.gene_start(idx) = G.gene_start(idx) - P.imputed_promoter_length;
  idx = grep('-',G.strand,1);
  G.gene_end(idx) = G.gene_end(idx) + P.imputed_promoter_length;
  G.footprint_including_promoter = G.gene_end-G.gene_start+1;
end

% load conservation if requested

if P.load_conservation29
  try
    [C Z] = load_track(['/cga/tcga-gsc/home/lawrence/db/' build '/conservation29']);
    G.tot_exon_len_cons29 = zeros(slength(G),1);
    for i=1:slength(G), for e=1:G.n_exons(i)
      if G.exon_ends{i}(e)>length(C{G.chr(i)})
        G.tot_exon_len_cons29(i) = nan;
        break;
      end
      G.tot_exon_len_cons29(i) = G.tot_exon_len_cons29(i) + ...
          sum(C{G.chr(i)}(G.exon_starts{i}(e):G.exon_ends{i}(e)));
    end,end
  catch me
    fprintf('Could not open conservation29 track for build %s\n',build);
  end
end

if P.load_conservation3
  try
    [C Z] = load_track(['/cga/tcga-gsc/home/lawrence/db/' build '/conservation3']);
    G.tot_exon_len_cons3 = zeros(slength(G),1);
    for i=1:slength(G), for e=1:G.n_exons(i)
      if G.exon_ends{i}(e)>length(C{G.chr(i)})
        G.tot_exon_len_cons3(i) = nan;
        break;
      end
      G.tot_exon_len_cons3(i) = G.tot_exon_len_cons3(i) + ...
          sum(C{G.chr(i)}(G.exon_starts{i}(e):G.exon_ends{i}(e)));
    end,end
  catch me
    fprintf('Could not open conservation3 track for build %s\n',build);
  end
end


