function process_hairpin_block_v7(blockno,mutfile,survdir,outdir,gcmin,gcmax)
% process_hairpin_block_v7(blockno,mutfile,survdir,outdir,gcmin,gcmax)
% --> adds gcmin, gcmax as parameters
% --> merges process_hairpin_block_vN and process_hairpin_block_vNa)

if ~exist('gcmin','var'), gcmin=0; end
if ~exist('gcmax','var'), gcmax=1; end
if ~isnumeric(gcmin), gcmin=str2double(gcmin); end
if ~isnumeric(gcmax), gcmax=str2double(gcmax); end
if gcmin<0 || gcmin>1 || gcmax<0 || gcmax>1, error('gcmin/gcmax should be 0-1'); end

if ~isnumeric(blockno), blockno=str2double(blockno); end
if blockno<1 || blockno>32, error('invalid blockno')'; end

load('hg19_genome_blocks.v1.0.mat','B');
chr=B.chr(blockno);st=B.st(blockno);en=B.en(blockno);

ede(outdir);
survfile  = [survdir  '/block' num2str(blockno) '.collapsed.mat'];
outfile = [outdir '/block' num2str(blockno) '.cts.mat'];
if exist(outfile,'file'), error('outfile already exists'); end; fprintf('BLOCK %d\n',blockno); save_textfile('in_process',outfile);

% LOAD MUTATION DATA
tic; load(mutfile,'X'); toc

% LOAD GENOMIC BLOCK from survey_hairpins_v*.m
tic; tmp=load(survfile);Z=tmp.X; toc; Z.site.chr = repmat(chr,slength(Z.site),1); Z=rmfield_if_exist(Z,'loop');
X.mut = reorder_struct(X.mut,X.mut.chr==chr & X.mut.pos>=st & X.mut.pos<=en);
tic; X.mut.zidx = listmap(X.mut.pos,Z.site.pos); toc

% ANNOTATE with replication timing, from hg19 windows reference file
%   (should have done this in survey_hairpins)
Z.site.reptime = nan(slength(Z.site),1);
load('hg19_windows.v1.0.mat','W'); W0=W;
if chr<24 % (don't have replication timing info for chrY
  W=reorder_struct(W.win,W.win.chr==chr);
  tic; Z.site.widx = mmw(Z.site,W); toc
  Z.site.reptime = nansub(W.rt_extra1,Z.site.widx);
end

% SITES TO USE

if gcmin>0 || gcmax<1
  % compute GC content in 100bp windows
  fprintf('gcmin = %d   gcmax = %d\n');
  tic; Z.site.gc100 = smooth((Z.site.ref==2 | Z.site.ref==3),100); toc
  Z.site.use = (Z.site.zone<3 & Z.site.gc100>=gcmin & Z.site.gc100<=gcmax);
else
  Z.site.use = (Z.site.zone<3);
end

% STEM STRENGTH = 3GC + 1AT
Z.site.stemstrength = Z.site.nbp + 2*Z.site.ngc;

% PATIENT SUBSETS (some are assoiated with a "patient-specific LEGO region")
subsets = {'apobec';'uv';'pole';'msi';'eso';'aging';'smoking';'brca';'cohort_nonapobec';'cohort_apobec';'cohort_A3A';'cohort_A3B'};
ns = length(subsets);

% TABULATE 
% make the following n and N tables:
%   COL  = looplen (3-11)
%   ROW  = stemstrength (0-26+)
%   PAGE = patient subsets (e.g. APOBEC, UV, POLE, MSI)

N = nan(27,9,ns); n = nan(27,9,ns);

% also measure dependence on replication time
Q=[];
Q.min=[100:100:1000]';
Q.max=[200:100:1000 inf]';
Q.N = nan(slength(Q),ns); Q.n=Q.N;

for si=1:ns,name=subsets{si};
  if strcmp(name,'apobec')  % Tp(C->G)
    puse = X.pat.frac_apobec>=0.5;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
    suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
  elseif strcmp(name,'uv')  % C->T
    puse = X.pat.frac_uv>=0.5 & X.pat.frac_apobec<0.1;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==4)|(X.mut.ref==3&X.mut.alt==1));
    suse = Z.site.use & (Z.site.ref==2 | Z.site.ref==3);
  elseif strcmp(name,'pole')  % all mutations
    puse = X.pat.frac_pole>=0.5 & X.pat.frac_apobec<0.1;
    muse = puse(X.mut.pat_idx);
    suse = Z.site.use;
  elseif strcmp(name,'msi')  % all mutations
    puse = X.pat.frac_msi>=0.5 & X.pat.frac_apobec<0.1;
    muse = puse(X.mut.pat_idx);
    suse = Z.site.use;
  elseif strcmp(name,'eso')  % A->C
    puse = X.pat.frac_eso>=0.5 & X.pat.frac_apobec<0.1;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==1&X.mut.alt==2)|(X.mut.ref==4&X.mut.alt==3));
    suse = Z.site.use & (Z.site.ref==1 | Z.site.ref==4);
  elseif strcmp(name,'aging')  % (C->T)pG
    puse = X.pat.frac_aging>=0.5 & X.pat.frac_apobec<0.1;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==4)|(X.mut.ref==3&X.mut.alt==1));
    suse = Z.site.use & ((Z.site.ref==2 & Z.site.plus1==3) | (Z.site.ref==3 & Z.site.minus0==2));
  elseif strcmp(name,'smoking')  % C->A
    puse = X.pat.frac_smoking>=0.5 & X.pat.frac_apobec<0.1;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==1)|(X.mut.ref==3&X.mut.alt==4));
    suse = Z.site.use & (Z.site.ref==2 | Z.site.ref==3);
  elseif strcmp(name,'brca')  % all mutations
    puse = X.pat.frac_brca>=0.5 & X.pat.frac_apobec<0.1;
    muse = puse(X.mut.pat_idx);
    suse = Z.site.use;
  elseif strcmp(name,'cohort_nonapobec')  % all mutations: need to repeat with Tp(C->G) only
    puse = X.pat.vanilla;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
    suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
  elseif strcmp(name,'cohort_apobec')  % Tp(C->G)
    puse = X.pat.apobec;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
    suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
  elseif strcmp(name,'cohort_A3A')  % Tp(C->G)
    puse = X.pat.apobec_most_a3a;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
    suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
  elseif strcmp(name,'cohort_A3B')  % Tp(C->G)
    puse = X.pat.apobec_most_a3b;
    muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
    suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
  else
    error('?');
  end

  fprintf('%20s %4d pat   %9d mut   %13d site\n',name,sum(puse),sum(muse));

  fld = [name '_ct'];
  Z.site.(fld) = histc(X.mut.zidx(muse),1:slength(Z.site));

  sidx0 = find(Z.site.use);
  for looplen=3:11
    sidx1 = sidx0(Z.site.looplen(sidx0)==looplen);
    for stemstrength=0:26
      if stemstrength<26
        sidx2 = sidx1(Z.site.stemstrength(sidx1)==stemstrength);
      else
        sidx2 = sidx1(Z.site.stemstrength(sidx1)>=stemstrength);
      end
      N(stemstrength+1,looplen-2,si) = sum(puse)*length(sidx2);
      n(stemstrength+1,looplen-2,si) = sum(Z.site.(fld)(sidx2));
    end
  end

  %% also measure dependence on replication time
  for i=1:slength(Q)
    sidx3 = sidx0(Z.site.reptime(sidx0)>=Q.min(i) & Z.site.reptime(sidx0)<Q.max(i));
    Q.N(i,si) = sum(puse)*length(sidx3);
    Q.n(i,si) = sum(Z.site.(fld)(sidx3));
  end

end % next si

save(outfile,'N','n','Q','subsets');

