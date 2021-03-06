%% CHANGE FOLLOWING LINE TO BE FULL PATH WHERE GITHUB REPO WAS COPIED TO %%

srcpath = '/data/rama/labMembers/al20/Lawrence-APOBEC/git/repo/';
addpath(genpath(srcpath));
maxNumCompThreads(8);
cd(srcpath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GENOMEWIDE SURVEY FOR HAIRPIN SITES
%
% --> Note, this code is slow, takes hours per chromosome.
% --> Instead of optimizing it, we submitted many jobs in parallel.
%
% It can either be run as a loop, or else "compiled" and run in parallel.

spliceflank=20;
outdir = 'v3a'; ede(outdir);
X=[]; tmp=load('hg19_genome_blocks.v1.0.mat','B'); X.block=tmp.B; X.gene = load_genes('hg19gencode',spliceflank);
X.loop=[];
X.loop.len = uint8([repmat(3,1,3) repmat(4,1,4) repmat(5,1,5) repmat(6,1,6) repmat(7,1,7) repmat(8,1,8) repmat(9,1,9) repmat(10,1,10) repmat(11,1,11)]');
X.loop.pos = uint8([1:3 1:4 1:5 1:6 1:7 1:8 1:9 1:10 1:11]');
X.loop.pos_flip = X.loop.len+1-X.loop.pos; X.loop.flip_idx = multimap(X.loop,X.loop,{'len','pos'},{'len','pos_flip'});
X.loop.tpos = (double(X.loop.pos)-1)-(double((X.loop.len+1))/2); X.loop.good = (abs(X.loop.tpos)<=0.5);
X.loop.pos_double = double(X.loop.pos); X.loop.len_double = double(X.loop.len); maxstem=250; window=maxstem+max(X.loop.len_double);
for bi=1:slength(X.block)
  outfile = [outdir '/block' num2str(bi) '.mat']; if exist(outfile,'file'), continue; end
  write_textfile('IN_PROGRESS',outfile);
  fprintf('BLOCK %d\n',bi);
  X.block.this = (1:slength(X.block))'==bi; chr=X.block.chr(bi); st=X.block.st(bi); en=X.block.en(bi); len=en-st+1;
  X.site=[]; X.site.pos=uint32((st:en)'); X.site.ref=uint8(listmap(upper(genome_region(chr,st,en,'hg19'))','ACGT'));
  X.site.gene = zeros(len,1,'uint16'); X.site.zone = zeros(len,1,'uint8'); % 0=IGR 1=intron 2=UTR 3=exon
  X.site.minus3  = [0;0;0;0;X.site.ref(1:end-4)];
  X.site.minus2  = [0;0;0;X.site.ref(1:end-3)];
  X.site.minus1  = [0;0;X.site.ref(1:end-2)];
  X.site.minus0  = [0;X.site.ref(1:end-1)];
  X.site.plus1   = [X.site.ref(2:end);0];
  X.site.plus2   = [X.site.ref(3:end);0;0];
  X.site.plus3   = [X.site.ref(4:end);0;0;0];
  X.site.plus4   = [X.site.ref(5:end);0;0;0;0];
  X.site.tpc = (X.site.ref==2 & X.site.minus0==4) | (X.site.ref==3 & X.site.plus1==1);
  for gi=1:slength(X.gene),if X.gene.chr(gi)~=chr, continue; end
    dst = max(1,X.gene.tx_start(gi)-st+1); den = min(len,X.gene.tx_end(gi)-st+1); X.site.gene(dst:den)=gi; X.site.zone(dst:den)=2;
    dst = max(1,X.gene.code_start(gi)-st+1); den = min(len,X.gene.code_end(gi)-st+1); X.site.gene(dst:den)=gi; X.site.zone(dst:den)=1;
  end
  for gi=1:slength(X.gene),if X.gene.chr(gi)~=chr, continue; end
    for ei=1:X.gene.n_exons(gi)
      dst = max(1,X.gene.exon_starts{gi}(ei)-st+1); den = min(len,X.gene.exon_ends{gi}(ei)-st+1); X.site.gene(dst:den)=gi; X.site.zone(dst:den)=3;
  end,end
  X.site.nbp = zeros(len,slength(X.loop),'uint8'); X.site.ngc = X.site.nbp;
  fprintf('Mb:');
  for i=window+1:len-window, if ~mod(i,1e6), fprintf(' %d/%d',i/1e6,ceil(len/1e6)); end
    for li=1:slength(X.loop), l = i-X.loop.pos_double(li); r = l+X.loop.len_double(li)+1; nbp=0; ngc=0;
      for k=1:maxstem, if X.site.ref(l)~=5-(X.site.ref(r)), break; end
        nbp=nbp+1; ngc=ngc+(X.site.ref(l)==2 || X.site.ref(l)==3); l=l-1; r=r+1;
      end
      X.site.nbp(i,li) = nbp; X.site.ngc(i,li) = ngc;
  end,end,fprintf('\n');
  ga = (X.site.ref==3 | X.site.ref==1); X.site.nbp(ga,:) = X.site.nbp(ga,X.loop.flip_idx); X.site.ngc(ga,:) = X.site.ngc(ga,X.loop.flip_idx); % flip info for ref=G/A
  tmp=X.site.plus1(ga); X.site.plus1(ga)=5-X.site.minus0(ga); X.site.minus0(ga)=5-tmp;
  tmp=X.site.plus2(ga); X.site.plus2(ga)=5-X.site.minus1(ga); X.site.minus1(ga)=5-tmp;
  tmp=X.site.plus3(ga); X.site.plus3(ga)=5-X.site.minus2(ga); X.site.minus2(ga)=5-tmp;
  tmp=X.site.plus4(ga); X.site.plus4(ga)=5-X.site.minus3(ga); X.site.minus3(ga)=5-tmp;
  save(outfile,'X','-v7.3');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optional: COMPILE AND RUN IN PARALEL

% COMPILE

% from bash shell:
%   cd /full/path/to/github/repo/directory
%   mkdir mcc
%   mkdir mcc/shp
%   mkdir mcc/shp/v7
%   cd mcc/shp/v7
%   cp /full/path/to/github/repo/directory/src/run_matlab.py .		% change to full path to repo directory
%   reuse .matlab-2016b                                             % change depending on environment; add matlab to environment path
%   mcc -m -C -I /full/path/to/github/repo/directory/ -I /full/path/to/github/repo/directory/ref \	% change to full path to repo directory
%       -d /full/path/to/github/repo/directory/mcc/shp/v7 survey_hairpins								% change to full path to current directory

% RUN IN PARALLEL

outdir = 'v3a'; ede(outdir);
mccdir = [srcdir '/mcc/shp/v7';
mccname = 'survey_hairpins';
load('hg19_genome_blocks.v1.0.mat','B');
cmds = {}; banner='survhp'; mem = '36g';
for blockno=1:slength(B)
  outfile = [outdir '/block' num2str(blockno) '.mat'];
  if exist(outfile,'file'), continue ;end
  cmds{end+1} = ['python run_matlab.py . ' mccname ' ' outdir ' ' num2str(blockno)];
end
cmdfile = [outdir '/cmd.txt']; save_lines(cmds,cmdfile); qcmd = ['qsubb ' cmdfile];
qcmd = [qcmd ' --pre="reuse .matlab-2016b;cd ' mccdir '"']; % change arguments to suit your cluster
qcmd = [qcmd ' -cwd -N ' banner ' -j y -l h_rt=48:00:00 -l h_vmem=' mem]; qcmd = [qcmd ' -o ' outdir];
qcmd = [qcmd char(10)]; qcmdfile = [outdir '/qcmd.txt']; save_textfile(qcmd,qcmdfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Continue here after running loop or parallel jobs:

% COLLAPSE AMBIGUOUS STRUCTURES
% collapse each position to its strongest stem (defined by 3GC+1AT), using shortest loop as tiebreaker
load('hg19_genome_blocks.v1.0.mat','B');
indir = 'v3a';
for bi=1:slength(B)
  infile=[indir '/block' num2str(bi) '.mat']; outfile = [indir '/block' num2str(bi) '.collapsed.mat'];
  if exist(outfile,'file'), continue; end
  fprintf('BLOCK %d ',bi);
  save_textfile('in_progress',outfile);
  load(infile,'X');
  X.site.dG3 = X.site.nbp+2*X.site.ngc; [tmp ord] = max(X.site.dG3,[],2);
  z = zeros(slength(X.site),1,'uint8'); X.site.best_nbp=z; X.site.best_ngc=z; X.site.best_looplen=z; X.site.best_looppos=z;
  for i=1:slength(X.loop), ii=(ord==i);
    X.site.best_nbp(ii)=X.site.nbp(ii,i); X.site.best_ngc(ii)=X.site.ngc(ii,i);
    X.site.best_looplen(ii)=X.loop.len(i); X.site.best_looppos(ii)=X.loop.pos(i);
  end
  X.site = rmfield(X.site,{'nbp','ngc','dG3'}); fs=grep('^best_',fieldnames(X.site));X.site=rename_fields(X.site,fs,regexprep(fs,'^best_',''));
  X = rmfield(X,'loop'); X.site=rmfield_if_exist(X.site,{'minus3','plus4'});
  save(outfile,'X','-v7.3');
end

% COMBINE TpCs FROM ALL BLOCKS
load('hg19_genome_blocks.v1.0.mat','B');
indir = 'v3a';
X=[];
for bi=1:slength(B),fprintf('BLOCK %d\n',bi);
  infile = [indir '/block' num2str(bi) '.collapsed.mat'];
  tmp=load(infile,'X');
  if bi==1, X.gene=tmp.X.gene; end
  X.site{bi} = reorder_struct(rmfield(tmp.X.site,{'tpc','minus0','plus3'}),tmp.X.site.tpc);
  X.site{bi}.chr = repmat(uint8(B.chr(bi)),slength(X.site{bi}),1);
end
X.site=concat_structs(X.site); X.site=order_field_first(X.site,'chr');
save([indir '/all.TpCs_only.mat'],'X','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WGS DATASET 
%

pat=makeapn(load_struct('FINAL_DATASET.MIN.v3.2.pat.txt'));
mut=arrayfun(@(x) makeapn(load_struct(['FINAL_DATASET.MIN.v3.2.mut.part' num2str(x) '.txt'])),[1:10]','uni',0);
mut=concat_structs(mut);
ttype=makeapn(load_struct('FINAL_DATASET.MIN.v3.2.ttype.txt'));

X=struct(); X.pat=pat; X.mut=mut; X.ttype=ttype;

save('FINAL_DATASET.MIN.v3.2.mat','X');

%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for reference:
%
% Sanger signatures 1-30

S = load_struct('sanger30.txt');
S = parsein(S,'channel','^(.)\[(.)>(.)\](.)$',{'l','f','t','r'});
bases = {'A','C','G','T'};S.f=listmap(S.f,bases); S.l=listmap(S.l,bases); S.r=listmap(S.r,bases); S.t=listmap(S.t,bases);
idx = find(S.f==4); S.f(idx)=5-S.f(idx); S.t(idx)=5-S.t(idx); tmp=5-S.r(idx); S.r(idx)=5-S.l(idx); S.l(idx)=tmp;
for i=1:slength(S),S.catname{i,1} = [bases{S.f(i)} ' in ' bases{S.l(i)} '_' bases{S.r(i)} ' ->' bases{S.t(i)}]; end
broad = generate_lego_96names(); S.ord = listmap(S.catname,broad); S = sort_struct(S,'ord')
X = {}; for i=1:30, X{i}=str2double(S.(['sig' num2str(i)]));end
save('sanger30.mat','X');
%%%

% show Lego plots of Sanger signatures 1-30

load('sanger30.mat','X');
legos(X,prefix(num2cellstr(1:30),'Signature '));
%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NMF ON FINAL COHORT
%
% k=8

load('FINAL_DATASET.MIN.v3.2.mat','X');

% make input for NMF, normalized to genomic territory
X.mut.patient = nansub(X.pat.name,X.mut.pat_idx);
P=[]; P.coding_only=false; X.nmf_input = make_input_for_nmf(X.mut,P);
X.mut = rmfield(X.mut,'patient');
X.nmf_input_rates = X.nmf_input;
load('breakdowns.v1.1.mat','K');
for i=1:slength(X.nmf_input_rates.chan)
  X.nmf_input_rates.chan.terr(i,1)=sum(K.Ng(K.f==X.nmf_input_rates.chan.f(i)&K.l==X.nmf_input_rates.chan.l(i)&K.r==X.nmf_input_rates.chan.r(i)));
end
X.nmf_input_rates.pat.nchan_counts = X.nmf_input_rates.pat.nchan;
X.nmf_input_rates.pat.nchan = bsxfun(@rdivide,X.nmf_input_rates.pat.nchan,X.nmf_input_rates.chan.terr');

k=8; randseed=1; X.nmf_rates = perform_nmf(X.nmf_input_rates,k,randseed);

figure(1);clf,display_nmf_legos(X.nmf_rates)                          % Manually inspect signatures...
names = {'APOBEC';'UV';'MSI';'smoking';'ESO';'aging';'POLE';'BRCA'};	% The FULL cohort had these 8 signatures.  Yours may vary.
ord = [1:8]                                                           % Note: POLE signature is absent from the LANZ cohort
X.nmf.factor.name = as_column(names);
X.nmf.chan = X.nmf_rates.chan;
X.nmf.chan.nmf = X.nmf_rates.chan.nmf(:,ord);
pord = listmap(X.pat.name,X.nmf_rates.pat.name);
X.pat.nchan_counts = X.nmf_rates.pat.nchan_counts(pord,:);
X.pat.nchan_rates = X.nmf_rates.pat.nchan(pord,:);
X.pat.nmf = X.nmf_rates.pat.nmf(pord,ord);
X.pat.nmf_norm = X.nmf_rates.pat.nmf_norm(pord,ord);
X=rmfield(X,{'nmf_input','nmf_input_rates','nmf_rates'});
X.pat=keep_fields(X.pat,{'name','cohort','ttype','ttype_long','ttype_idx','nmut','nchan_counts','nchan_rates','nmf','nmf_norm'});
X.pat.frac_apobec = X.pat.nmf_norm(:,1);
X.pat.frac_uv = X.pat.nmf_norm(:,2);
X.pat.frac_msi = X.pat.nmf_norm(:,3);
X.pat.frac_smoking = X.pat.nmf_norm(:,4);
X.pat.frac_eso = X.pat.nmf_norm(:,5);
X.pat.frac_aging = X.pat.nmf_norm(:,6);
X.pat.frac_pole = X.pat.nmf_norm(:,7);
X.pat.frac_brca = X.pat.nmf_norm(:,8);

figure(1);clf,display_nmf_legos(X.nmf)

save('FINAL_DATASET.v2.1.mat','X');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEGO PLOTS

load('FINAL_DATASET.v2.1.mat','X');

figure(1)
ede('finalfigs/sampledata/v1');
P=[];P.catnames=X.nmf.chan.name;P.add_sanger_plot=true;P.use_lighter_blue=true;
for i=1:slength(X.nmf.factor),disp(X.nmf.factor.name{i});
  clf,lego(X.nmf.chan.nmf(:,i),P);
  set(gcf,'papersize',[3 3],'paperposition',[0.2 0.2 2.6 2.6]);
  print_to_file(['finalfigs/sampledata/v1/lego.' X.nmf.factor.name{i} '.pdf']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add more info to mutation object
% --> define patient cohorts

load('FINAL_DATASET.v2.1.mat','X');
load('pats_with_ABdel_classifications.mat','pat');
X.pat = parsein(X.pat,'name','^[^:]+:(.+)$','name0');
X.pat.ABdel = mapacross(X.pat.name0,pat.name,pat.isdel);
X.pat.nmut = histc(X.mut.pat_idx,1:slength(X.pat));
X.pat=order_fields_first(X.pat,{'name','name0','cohort','nmut'});
tic;X.mut = add_llftrr(X.mut,[],[srcpath 'ref/pentamer/']);toc; % ~3min
X.mut = rmfield(X.mut,{'name','num'});
c2gt = (X.mut.f==2 & X.mut.t~=1);
X.pat.nC = histc(X.mut.pat_idx(c2gt),1:slength(X.pat));
tca2gt = (c2gt & X.mut.l==4 & X.mut.r==1);
X.pat.nRTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==1 | X.mut.ll==3)),1:slength(X.pat));
X.pat.nYTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==2 | X.mut.ll==4)),1:slength(X.pat));
X.pat.RTCA_C = X.pat.nRTCA./X.pat.nC; X.pat.YTCA_C = X.pat.nYTCA./X.pat.nC;
x = X.pat.RTCA_C; y = X.pat.YTCA_C;
X.pat.msupe_neg = (X.pat.frac_msi<0.1&X.pat.frac_smoking<0.1&X.pat.frac_uv<0.1&X.pat.frac_pole<0.1&X.pat.frac_eso<0.1);
X.pat.vanilla = (X.pat.frac_apobec<0.02 & X.pat.msupe_neg); 
X.pat.apobec =  (X.pat.frac_apobec>=0.1 & X.pat.msupe_neg);   
X.pat.apobec_most_a3a = X.pat.apobec & y>0.281;  
X.pat.apobec_most_a3b = X.pat.apobec & x>0.05 & y<0.116 & (x./(y+0.0315)>=0.655);  
save('FINAL_DATASET.v2.2.mat','X','-v7.3');
%%%%%%%%%%%%%%%


% output table of patients
load('FINAL_DATASET.v2.2.mat','X');
pat = X.pat;
pat = sort_struct(pat,{'cohort','ttype','name0'});
pat = keep_fields(pat,{'cohort','name0','ttype','ttype_long','nmut','ABdel','frac_apobec','frac_uv','frac_pole','frac_msi','frac_smoking','frac_eso','frac_aging','frac_brca',...
                    'msupe_neg','vanilla','apobec','RTCA_C','YTCA_C','apobec_most_a3a','apobec_most_a3b'});
pat=rename_field(pat,'name0','name');
ede('supp_data_files');
save_struct(pat,'supp_data_files/patient_list.txt');
%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAP final dataset onto genomic windows

load('FINAL_DATASET.v2.2.mat','X');
load('hg19_windows.v1.0.mat','W');
X.gene=W.gene; X.win=W.win; clear W
X.mut.zone_idx = get_context(X.mut.chr,X.mut.pos,[srcpath 'ref/zone/']);
X.mut.win_idx = mmw(X.mut,X.win);
puse = (1:slength(X.pat)); muse = find(ismember(X.mut.pat_idx,puse) & X.mut.zone_idx~=3); % nonexon only
X.win.nmut_nonexon_newdata = histc(X.mut.win_idx(muse),1:slength(X.win));
X.win.mutrate_nonexon_newdata = X.win.nmut_nonexon_newdata./(X.win.cov_nonexon*slength(X.pat));
X.win.mutrate_nonexon_newdata(isinf(X.win.mutrate_nonexon_newdata))=nan;  % (smooth behaves differently on nan vs. inf)
save('FINAL_DATASET_on_windows.v2.2.mat','X','-v7.3');
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIGURE 1a
% large-scale covariates & mutation rate

load('FINAL_DATASET_on_windows.v2.2.mat','X');

chrlen = load_chrlen('hg19'); cen = load_cen('hg19'); window=100000;
chr=17; st=1; en=81e6; ptel=1e5; pcen=1e6; qcen=1e5; qtel=1e6;
Wc = reorder_struct(X.win,X.win.chr==chr);
dat1 = 150+smooth(Wc.rt_extra1,20);
dat2 = 150+max(200,1100-smooth(Wc.expr/500,16));
dat3 = -200+smooth(2.6e8*Wc.mutrate_nonexon_newdata,20);
dat4 = 600-smooth(Wc.comp*11500,20);
datz = [dat1 dat2 dat3 dat4];

figure(1);clf;hold on;fontsize=26;
left = max(1,min(chrlen(chr),window*round(st/window))); right= max(1,min(chrlen(chr),window*round(en/window)));
bins = ceil((left:window:right)/window); pos=Wc.pos(bins);
cenidx = find(left:window:right>=cen(chr,1)-pcen&qcen+cen(chr,2)>=left:window:right); % cen mask
telidx = find(ptel>left:window:right|left:window:right>chrlen(chr)-qtel); % tel mask
bad = union(cenidx,telidx); datz(bad,:) = nan;
plot(pos/1e6,datz(bins,1),'-','linewidth',6,'color',[0 0 1]);         % reptime
plot(pos/1e6,datz(bins,2),'-','linewidth',6,'color',[0.2 0.7 0.2]);   % expression
plot(pos/1e6,datz(bins,3),'-','linewidth',6,'color',[1 0 0]);         % mutrate nonexon
plot(pos/1e6,datz(bins,4),'-','linewidth',6,'color',[1 0.75 0]);      % HiC
xlim([-1+pos(1)/1e6 pos(end)/1e6+1]); ylim([-750 2400]); ff; set(gca,'visible','off');
genes = {'TP53';'BRCA1';'ERBB2';'NF1';'KIF2B';'ABCA5';'DNAH9';'ASIC2'};
gbary = [120;    120;    120;    120;   1790;   1450;  1600;   1330];
gbarh = [ 30;     30;     30;     30;     30;     30;    30;     30];
gtxty = [-45;    -45;    -45;    -45;   1980;   1650;  1800;   1520];
gidx = listmap(genes,X.gene.name); x=[];y=[];txt={};
for i=1:length(gidx),gi=gidx(i);
  if X.gene.chr(gi)~=chr, continue; end
  wst = X.gene.tx_start(gi)/1e6; wen = X.gene.tx_end(gi)/1e6;
  rectangle('position',[wst gbary(i) wen-wst gbarh(i)],'edgecolor',[0 0 0],'facecolor',[0 0 0],'linewidth',6);
  x(end+1,1)=mean([wst wen]); y(end+1,1)=gtxty(i); txt{end+1,1}=genes{i};
end, x(3)=x(3)-2; x(2)=x(2)+2; text(x,y,txt,'fontsize',16,'hor','cen');

set(gcf,'papersize',[13 5],'paperposition',[0 0.05 12.9 4.8]);
print_to_file('finalfigs/sampledata/v1/fig1a.pdf');
%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A3A vs. A3B PLOTS

load('FINAL_DATASET.v2.2.mat','X');

% Fig. S4a 
figure(3);clf,hold on
x = X.pat.RTCA_C; y = X.pat.YTCA_C;
clr = repmat([0.8 0.8 0.8],slength(X.pat),1); sz=10*ones(slength(X.pat),1);
van = (X.pat.vanilla); clr(van,:)=repmat([0.2 0.7 0.2],sum(van),1);sz(van)=30;
apo = (X.pat.apobec); clr(apo,:)=repmat([0 0 0],sum(apo),1);sz(apo)=30;
a3b = (X.pat.apobec_most_a3b);clr(a3b,:)=repmat([0 0 1],sum(a3b),1); sz(a3b)=50;
a3a = (X.pat.apobec_most_a3a);clr(a3a,:)=repmat([1 0 0],sum(a3a),1); sz(a3a)=50;
abd = (X.pat.ABdel==1);
scatter(x,y,sz,clr,'filled');
scatter(x(van),y(van),sz(van),clr(van,:),'filled');
scatter(x(apo),y(apo),sz(apo),clr(apo,:),'filled');
scatter(x(a3b),y(a3b),sz(a3b),clr(a3b,:),'filled');
scatter(x(a3a),y(a3a),sz(a3a),clr(a3a,:),'filled');
scatter(x(abd),y(abd),70,[1 0.7 0],'filled');scatter(x(abd),y(abd),90,[0 0 0]);scatter(x(abd),y(abd),sz(abd),clr(abd,:),'filled');
ff;set(gca,'fontsize',20,'ytick',0:0.1:4);xlabel('APOBEC3B character','fontsize',30); ylabel('APOBEC3A character','fontsize',30);xlim([0 0.14]);
set(gcf,'papersize',[9 8],'paperposition',[0.2 0.2 9-0.4 8-0.4]);
print_to_file('finalfigs/sampledata/v1/fig_S4a.pdf');

% Fig. S4c
ede('finalfigs/sampledata/v1');
mut_nonapobec = reorder_struct(X.mut,X.pat.vanilla(X.mut.pat_idx));
mut_apobec = reorder_struct(X.mut,X.pat.apobec(X.mut.pat_idx));
mut_apobec_most_a3b = reorder_struct(X.mut,X.pat.apobec_most_a3b(X.mut.pat_idx));
mut_apobec_most_a3a = reorder_struct(X.mut,X.pat.apobec_most_a3a(X.mut.pat_idx));
figure(1),clf,set(gcf,'papersize',[3 3],'paperposition',[0.2 0.2 2.6 2.6]);
P=[];P.display='rates';P.coverage='genome';P.add_sanger_plot=true;P.use_lighter_blue=true;
clf,legomaf(mut_nonapobec,P);print_to_file('finalfigs/sampledata/v1/lego.cohort.non_APOBEC.pdf');
clf,legomaf(mut_apobec,P);print_to_file('finalfigs/sampledata/v1/lego.cohort.APOBEC.pdf');
clf,legomaf(mut_apobec_most_a3b,P);print_to_file('finalfigs/sampledata/v1/lego.cohort.A3B.pdf');
clf,legomaf(mut_apobec_most_a3a,P);print_to_file('finalfigs/sampledata/v1/lego.cohort.A3A.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
%
% PROCESS EACH GENOMIC BLOCK to calculate large- and mesoscale response curves
%
% --> This code is slow.
% --> Instead of optimizing it, we submitted many jobs in parallel.
%
% It can either be run as a loop, as below, or else "compiled" and run in parallel.
%

mutfile = 'FINAL_DATASET.v2.2.mat';
survdir = [srcdir 'wgs/v3/run12/'];
survfile  = [survdir  '/block' num2str(blockno) '.collapsed.mat'];
gcmin=0.4; gcmax=0.6;
outdir = 'wgs/v3/run12'; ede(outdir);

% LOAD MUTATION DATA
tic; load(mutfile,'X'); toc
Xo=X

load('hg19_genome_blocks.v1.0.mat','B');      

for blockno=1:32

  chr=B.chr(blockno);st=B.st(blockno);en=B.en(blockno);

  outfile = [outdir '/block' num2str(blockno) '.cts.mat'];
  if exist(outfile,'file'), error('outfile already exists'); end; fprintf('BLOCK %d\n',blockno); save_textfile('in_process',outfile);

  X=Xo

  % LOAD GENOMIC BLOCK from survey_hairpins_v*.m
  tic; tmp=load(survfile);Z=tmp.X; toc; Z.site.chr = repmat(chr,slength(Z.site),1); Z=rmfield_if_exist(Z,'loop');
  X.mut = reorder_struct(X.mut,X.mut.chr==chr & X.mut.pos>=st & X.mut.pos<=en);
  tic; X.mut.zidx = listmap(X.mut.pos,Z.site.pos); toc

  % ANNOTATE with replication timing, from hg19 windows reference file
  Z.site.reptime = nan(slength(Z.site),1);
  load('hg19_windows.v1.0.mat','W'); W0=W;
  if chr<24 % (don't have replication timing info for chrY
    W=reorder_struct(W.win,W.win.chr==chr);
    W.chr = repmat(chr,slength(W),1);
    W.st=W.pos+1; W.en=W.pos+window;
    tic; Z.site.widx = mmw(Z.site,W); toc
    Z.site.reptime = nansub(W.rt_extra1,Z.site.widx);
  end

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

  % PATIENT SUBSETS 
  subsets = {'apobec';'uv';'pole';'msi';'eso';'aging';'smoking';'brca';'cohort_nonapobec';'cohort_apobec';'cohort_A3A';'cohort_A3B'};
  ns = length(subsets);

  % TABULATE
  % make the following n and N tables:
  %   COL  = looplen (3-11)
  %   ROW  = stemstrength (0-26+)
  %   PAGE = patient subsets (APOBEC, UV, POLE, MSI, etc.)

  N = nan(27,9,ns); n = nan(27,9,ns);

  % also measure dependence on replication time
  Q=[]; % reptime deciles
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
    elseif strcmp(name,'cohort_nonapobec') 
      puse = X.pat.vanilla;
      muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
      suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
    elseif strcmp(name,'cohort_apobec')
      puse = X.pat.apobec;
      muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
      suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
    elseif strcmp(name,'cohort_A3A') 
      puse = X.pat.apobec_most_a3a;
      muse = puse(X.mut.pat_idx) & ((X.mut.ref==2&X.mut.alt==3)|(X.mut.ref==3&X.mut.alt==2));
      suse = Z.site.use & ((Z.site.ref==2 & Z.site.minus0==4) | (Z.site.ref==3 & Z.site.plus1==1));
    elseif strcmp(name,'cohort_A3B') 
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

end % next blockno


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optional: COMPILE AND RUN IN PARA

% COMPILE

% from bash shell:
%   cd /full/path/to/github/repo/directory
%   mkdir mcc
%   mkdir mcc/phb
%   mkdir mcc/phb/v7
%   cd mcc/phb/v7
%   cp /full/path/to/github/repo/directory/src/run_matlab.py .		% change to full path to repo directory
%   echo "-Xmx16g" > java.opts
%   reuse .matlab-2016b                                             % change depending on environment; add matlab to environment path
%   mcc -m -C -I /full/path/to/github/repo/directory/ -I /full/path/to/github/repo/directory/ref \	% change to full path to repo directory
%       -d /full/path/to/github/repo/directory/mcc/phb/v7 process_hairpin_block_v7					% change to full path to current directory

% RUN IN PARALLEL

mutfile = 'FINAL_DATASET.v2.2.mat';
survdir = [srcdir 'wgs/v3/run12/'];
gcmin=0.4; gcmax=0.6;
outdir = 'wgs/v3/run12'; ede(outdir);
mccdir = [srcdir '/mcc/phb/v7';
mccname = 'process_hairpin_block_v7';
cmds={}; banner='phb';
for blockno=1:32
  cmds{end+1} = ['python run_matlab.py . ' mccname ' ' num2str(blockno) ' ' mutfile ' ' survdir ' ' outdir ' ' num2str(gcmin) ' ' num2str(gcmax)];
end
cmdfile = [outdir '/cmd.txt']; save_lines(cmds,cmdfile); qcmd = ['qsubb ' cmdfile];
qcmd = [qcmd ' --pre="reuse .matlab-2016b;cd ' mccdir '"']; 	% change arguments to suit your cluster
qcmd = [qcmd ' -cwd -N ' banner ' -j y -l h_rt=2:00:00 -l h_vmem=60g']; qcmd = [qcmd ' -o ' outdir]; 
qcmd = [qcmd char(10)]; qcmdfile = [outdir '/qcmd.txt']; save_textfile(qcmd,qcmdfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Continue here after running loop or parallel jobs:

% gather
n=0;N=0;Q=[];nn=[];NN=[]; sss=[0:26]';
indir = 'wgs/v3/run12';
for blockno=1:32, infile = [indir '/block' num2str(blockno) '.cts.mat']; tmp=load(infile,'N','n','Q','subsets');
  n=n+tmp.n; N=N+tmp.N; if blockno==1, subsets=tmp.subsets; Q=tmp.Q; else Q.N=Q.N+tmp.Q.N; Q.n=Q.n+tmp.Q.n;
end,end
ssmin = [0 4:2:16 20]; ssmax = [ssmin(2:end)-1 inf];
for i=1:length(subsets), row=0;
  for ssi=1:length(ssmin), row=row+1;
    j=(sss>=ssmin(ssi) & sss<=ssmax(ssi));
    nn(row,:,i) = sum(n(j,:,i),1); NN(row,:,i) = sum(N(j,:,i),1);
end,end, r = nn./NN;
nnn=squeeze(sum(nn(:,1:3,:),2)); NNN=squeeze(sum(NN(:,1:3,:),2)); [r sd] = ratio_and_sd(nnn,NNN);
rmd = nanmedian(r(1:2,:),1); r = bsxfun(@rdivide,r,rmd); sd = bsxfun(@rdivide,sd,rmd);
save('wgs/v3/run12/stats_all.v1.1.mat','nn','NN','Q','subsets');

% print numbers
prn(ssmin',ssmax',r)

%%%

load('stats_all.v1.1.mat','NN','nn','Q','subsets');

% hairpin-potential decile plots
figure(6),clf
nnn=squeeze(sum(nn(:,1:3,:),2)); NNN=squeeze(sum(NN(:,1:3,:),2)); [r sd] = ratio_and_sd(nnn,NNN);
rmd = nanmedian(r(1:2,:),1); r = bsxfun(@rdivide,r,rmd); sd = bsxfun(@rdivide,sd,rmd);
for i=1:length(subsets)
  subplot(ceil(length(subsets)/2),2,i)
  b = barweb(r(:,i),r(:,i)-1.96*sd(:,i),r(:,i)+1.96*sd(:,i));colormap(jet);ff
  for j=1:length(b.bars),set(b.bars(j),'linewidth',2);end
  ylim([0 11]);set(gca,'ytick',0:2:10);
  set(gca,'linewidth',2,'fontsize',26,'xtick',[]);
  line(xlim,[1 1],'color',[0 0 0],'linestyle','--');
  title(regexprep(subsets{i},'_','\\_'),'fontsize',30);
  if i==5, text(0.25,5.15,'relative mutation rate','rot',90,'fontsize',30); end
end
% replication-timing decile plots
figure(7),clf
[r sd] = ratio_and_sd(Q.n,Q.N);
med = nanmedian(r,1);r=bsxfun(@rdivide,r,med);sd=bsxfun(@rdivide,sd,med);
for i=1:length(subsets)
  subplot(ceil(length(subsets)/2),2,i)
  b = barweb(r(:,i),r(:,i)-1.96*sd(:,i),r(:,i)+1.96*sd(:,i));colormap(parula);ff
  for j=1:length(b.bars),set(b.bars(j),'linewidth',2);end
  set(gca,'linewidth',2,'fontsize',26,'xtick',[]);
  line(xlim,[1 1],'color',[0 0 0],'linestyle','--');
  title(subsets{i},'fontsize',30,'interp','none');
end

%(vector-graphics): hairpin-potential decile plots
load('stats_all.v1.1.mat','NN','nn','Q','subsets');
figure(2)
nnn=squeeze(sum(nn(:,1:3,:),2)); NNN=squeeze(sum(NN(:,1:3,:),2)); [r sd] = ratio_and_sd(nnn,NNN);
rmd = nanmedian(r(1:2,:),1); r = bsxfun(@rdivide,r,rmd); sd = bsxfun(@rdivide,sd,rmd);
effs=[]; corrs=[]; pvals=[];
for i=1:length(subsets)
  clf, b = barweb(r(:,i),r(:,i)-1.96*sd(:,i),r(:,i)+1.96*sd(:,i));colormap(jet);
  effs(i,1) = r(end,i)/r(1,i); [corrs(i,1) pvals(i,1)] = corrr((1:size(r,1))',r(:,i));
  for j=1:length(b.bars),set(b.bars(j),'linewidth',1);end;xlim([0.55 1.45]);ylim([0 11]);
  set(gca,'linewidth',1.5,'fontsize',16,'xtick',[]);yl=ylim;
  line(xlim,[1 1],'color',[0 0 0],'linestyle',':');
  set(gca,'ytick',[0:2:10],'ticklength',[0.015 0.015]);
  set(gcf,'papersize',[4.6 4.2],'paperposition',[0.4 0.4 3.8 3.3],'color',[1 1 1]);
  print_to_file(['finalfigs/sampledata/v1/hairpins.' subsets{i} '.pdf'],600);
end
pr(subsets,effs,corrs,pvals);

%(vector-graphics): replication-timing decile plots
load('wgs/v3/run11/stats_all.v1.1.mat','NN','nn','Q','subsets');
figure(2)
[r sd] = ratio_and_sd(Q.n,Q.N);
med = nanmedian(r,1); r=bsxfun(@rdivide,r,med);sd=bsxfun(@rdivide,sd,med);
effs=[]; corrs=[]; pvals=[];
for i=1:length(subsets)
  clf, b = barweb(r(:,i),r(:,i)-1.96*sd(:,i),r(:,i)+1.96*sd(:,i));colormap(parula);
  effs(i,1) = r(end,i)/r(1,i); [corrs(i,1) pvals(i,1)] = corrr((1:size(r,1))',r(:,i));
  for j=1:length(b.bars),set(b.bars(j),'linewidth',1);end;xlim([0.55 1.45]); if i<=4, ylim([0 1.5]); end
  set(gca,'linewidth',1.5,'fontsize',16,'xtick',[]);yl=ylim;if yl(2)<1.5, ylim([0 1.5]); end
  line(xlim,[1 1],'color',[0 0 0],'linestyle',':');
  set(gca,'ytick',[0:0.5:5],'ticklength',[0.015 0.015]);
  set(gcf,'papersize',[4.6 4.2],'paperposition',[0.4 0.4 3.8 3.3],'color',[1 1 1]);
  print_to_file(['finalfigs/sampledata/v1/reptime.' subsets{i} '.pdf'],600);
end
pr(subsets,log2(effs),corrs,pvals);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAP MUTATIONS TO LIST OF TpCs

tpcfile = 'all.TpCs_only.mat';
mutfile = 'FINAL_DATASET.v2.2.mat';
tic;load(tpcfile,'X');toc % 40 sec
tic;tmp=load(mutfile);toc % 10 sec
fs={'pat','ttype','nmf','mut'};for fi=1:length(fs),f=fs{fi};X.(f)=tmp.X.(f);end
%tic;X.mut.sidx = multimap(X.mut,X.site,{'chr','pos'});toc % 12 min
[lia,X.mut.sidx]=ismember(table([X.mut.chr X.mut.pos]),table([X.site.chr X.site.pos]),'rows');
X.mut.sidx(X.mut.sidx==0)=NaN;
%fld1=fieldnames(X.mut); fld2=setdiff(fieldnames(X.site),fld1);
%for j=1:length(fld2) X.mut.(fld2{j})=repmat(NaN,slength(X.mut),1); X.mut.(fld2{j})(lia)=X.site.(fld2{j})(locb(lia)); end
save('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figs.3b,S6

% TABULATE 

tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec to load --> 9.3Gb in mem
puse = (X.pat.frac_apobec>=0.5); muse = puse(X.mut.pat_idx); S=X.site; S.ct = histc(X.mut.sidx(muse),1:slength(S)); S.ss = S.nbp+2*S.ngc;
keep = (S.zone<3); S = reorder_struct(keep_fields(S,{'ct';'looplen';'looppos';'ss'}),keep); ii = (S.ss<5); nohprate = mean(S.ct(ii));
bin=[];bin.min=[0 7:2:17 19 18]';bin.max=[6 8:2:18 inf inf]';
bin.label = num2cellstr(bin.min); bin.label{end-1}(end+1)='+'; bin.label{end}(end+1)='+'; bin.label{1}=['<' bin.label{2}]; Qs={};
for looplen=3:11,disp(looplen);S1 = reorder_struct(S,S.looplen==looplen); Q=bin;
  for looppos=1:looplen, S2 = reorder_struct(S1,S1.looppos==looppos);
    for i=1:slength(Q), ii = (S2.ss>=Q.min(i) & S2.ss<Q.max(i));
      Q.N(i,looppos) = sum(ii); Q.n(i,looppos) = sum(S2.ct(ii));
  end,end, [Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N);Q.relrate = Q.rate/nohprate; Q.relsd = Q.sd/nohprate; Qs{looplen-2}=Q;
end % ~2 min
save('FINAL.Qs.v3a.1.mat','Qs');

% PLOT for Fig. 3b
load('FINAL.Qs.v3a.1.mat','Qs');
figure(4),clf,looplen=4;Q=Qs{looplen-2};Q=reorder_struct(Q,1:slength(Q)-1);rate=Q.relrate;sd=Q.relsd;ratehi=rate+1.96*sd;ratelo=rate-1.96*sd;
b=barweb(rate,ratelo,ratehi);ff;colormap(parula);set(gca,'xtick',1:slength(Q),'xticklabel',Q.label,'fontsize',14,'linewidth',1.5);
for bi=1:length(b.bars),set(b.bars(bi),'linewidth',1);end;xlim([0.3 slength(Q)+0.7]);line(xlim,[1 1],'linestyle',':','color',[0 0 0]);
set(gcf,'papersize',[10 7.5],'paperposition',[0.2 0.2 10-0.4 7.5-0.4]);
print_to_file('finalfigs/sampledata/v1/fig.3b.pdf');

% PLOT for Fig. S6
load('FINAL.Qs.v3a.1.mat','Qs');
figure(5),clf,looplen=3;Q=Qs{looplen-2};Q=reorder_struct(Q,1:slength(Q)-1);rate=Q.relrate;sd=Q.relsd;ratehi=rate+1.96*sd;ratelo=rate-1.96*sd;
b=barweb(rate,ratelo,ratehi);ff;colormap(parula);set(gca,'xtick',1:slength(Q),'xticklabel',Q.label,'fontsize',14,'linewidth',1.5);
for bi=1:length(b.bars),set(b.bars(bi),'linewidth',1);end;xlim([0.3 slength(Q)+0.7]);line(xlim,[1 1],'linestyle',':','color',[0 0 0]);
set(gcf,'papersize',[10 7.5],'paperposition',[0.2 0.2 10-0.4 7.5-0.4]);
print_to_file('finalfigs/sampledata/v1/fig.S6.pdf');

% PLOT for Fig. 3d
load('FINAL.Qs.v3a.1.mat','Qs');
figure(6),clf,hold on,cmap=parula;maxlooplen=8;x=0;xt=[];
for looplen=3:maxlooplen,Q=Qs{looplen-2};Q=reorder_struct(Q,slength(Q));x=x+2; xt(end+1)=x+0.5+looplen/2;
  rate=Q.relrate;sd=Q.relsd;ratelo=rate-1.96*sd;ratehi=rate+1.96*sd;
  for looppos=1:looplen,y=rate(looppos);ylo=ratelo(looppos);yhi=ratehi(looppos);
    x=x+1; bar(x,y,1,'facecolor',cmap(1+floor((size(cmap,1)-1)*(looppos-1)/(looplen-1)),:),'linewidth',1);
    ebw=0.25;line([x x],[ylo yhi],'color',[0 0 0],'linewidth',0.75);line(x+[-1 1]*ebw,[1;1]*[ylo yhi],'color',[0 0 0],'linewidth',1);
end,end
xlim([0.46 x+1.9]);ylim([0 49]);set(gca,'ytick',0:10:40,'xtick',xt,'xticklabel',3:maxlooplen,'fontsize',14,'linewidth',1.5);ff
line(xlim,[1 1],'linestyle',':','color',[0 0 0]);colorbar('position',[0.6 0.8 0.3 0.02],'orientation','horizontal');
set(gcf,'papersize',[9.5 7.5],'paperposition',[0.2 0.2 9.5-0.4 7.5-0.4]);
print_to_file('finalfigs/sampledata/v1/fig.3d.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figs.3f,S8

% TABULATE 

tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec to load --> 9.3Gb in mem
puse = (X.pat.frac_apobec>=0.5); muse = puse(X.mut.pat_idx);
So=X.site; So.ct = histc(X.mut.sidx(muse),1:slength(So)); So.ss = So.nbp+2*So.ngc; ii = (So.ss<5); nohprate = mean(So.ct(ii));
So=reorder_struct(keep_fields(So,{'looplen';'looppos';'ss';'minus2';'minus1';'plus1';'plus2';'ct'}),...
  So.zone<3 & ((So.looplen==3&So.looppos==3)|(So.looplen==4&So.looppos==3)|(So.looplen==4&So.looppos==4)));

bbb={'A';'C';'G';'T'}; rc='';rc('ACGT')='TGCA';
S=reorder_struct(So,So.looplen==3 & So.looppos==3); divs=[0 5:2:19];  % triloops: required dG3>=19
context=[];Q=[]; Q.min = as_col(divs); Q.max = [Q.min(2:end);inf]; Q.label = num2cellstr(Q.min); Q.label{end}(end+1)='+';
j=0;,for plus1=1:4,S1=reorder_struct(S,S.plus1==plus1); for minus1=1:4, S2=reorder_struct(S1,S1.minus1==minus1);
    j=j+1; context.name{j,1} = [bbb{minus1} 'TC' bbb{plus1}];
    for i=1:slength(Q), ii = (S2.ss>=Q.min(i) & S2.ss<Q.max(i));
      Q.N(i,j) = sum(ii); Q.n(i,j) = sum(S2.ct(ii));
end,end,end
[Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N);Q.relrate = Q.rate/nohprate; Q.relsd = Q.sd/nohprate; Q.cv=Q.relsd./Q.relrate;
QL3=Q;contextL3=context;
S=reorder_struct(So,So.looplen==4 & So.looppos==4); divs=[0 5:2:13];  % tetraloops: required dG3>=13
context=[];Q=[]; Q.min = as_col(divs); Q.max = [Q.min(2:end);inf]; Q.label = num2cellstr(Q.min); Q.label{end}(end+1)='+';
j=0;for plus1=1:4,S1=reorder_struct(S,S.plus1==plus1); for minus1=1:4, S2=reorder_struct(S1,S1.minus1==minus1);
    for minus2=1:4, S3=reorder_struct(S2,S2.minus2==minus2);
      j=j+1; context.name{j,1} = [bbb{minus2} bbb{minus1} 'TC' bbb{plus1}];
      for i=1:slength(Q), ii = (S3.ss>=Q.min(i) & S3.ss<Q.max(i));
        Q.N(i,j) = sum(ii); Q.n(i,j) = sum(S3.ct(ii));
end,end,end,end
[Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N);Q.relrate = Q.rate/nohprate; Q.relsd = Q.sd/nohprate; Q.cv=Q.relsd./Q.relrate;
Q4=Q;context4=context;
S=reorder_struct(So,So.looplen==4 & So.looppos==3); divs=[0 5:2:13];  % tetraloops: required dG3>=13
context=[];Q=[]; Q.min = as_col(divs); Q.max = [Q.min(2:end);inf]; Q.label = num2cellstr(Q.min); Q.label{end}(end+1)='+';
j=0;for plus1=1:4,S1=reorder_struct(S,S.plus1==plus1); for minus1=1:4, S2=reorder_struct(S1,S1.minus1==minus1);
    for minus2=1:4, S3=reorder_struct(S2,S2.minus2==minus2);
      j=j+1; context.name{j,1} = [bbb{minus2} bbb{minus1} 'TC' bbb{plus1}];
      for i=1:slength(Q), ii = (S3.ss>=Q.min(i) & S3.ss<Q.max(i));
        Q.N(i,j) = sum(ii); Q.n(i,j) = sum(S3.ct(ii));
end,end,end,end
[Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N);Q.relrate = Q.rate/nohprate; Q.relsd = Q.sd/nohprate; Q.cv=Q.relsd./Q.relrate;
Q3=Q;context3=context;
Z3=reorder_struct(Q3,6);Z3=rmfield(Z3,{'min','max','label'});fs=fieldnames(Z3);for fi=1:length(fs),f=fs{fi};Z3.(f)=Z3.(f)';end
Z4=reorder_struct(Q4,6);Z4=rmfield(Z4,{'min','max','label'});fs=fieldnames(Z4);for fi=1:length(fs),f=fs{fi};Z4.(f)=Z4.(f)';end
Z4.looplen=repmat(4,slength(Z4),1);Z4.looppos=repmat(4,slength(Z4),1);Z4.context=context4.name;
Z3.looplen=repmat(4,slength(Z3),1);Z3.looppos=repmat(3,slength(Z3),1);Z3.context=context4.name;
for i=1:slength(Z3),Z3.context{i}=[Z3.context{i}(1) '(' Z3.context{i}(2:5) ')' rc(Z3.context{i}(1))];end
for i=1:slength(Z4),Z4.context{i}=[rc(Z4.context{i}(5)) '(' Z4.context{i}(1:4) ')' Z4.context{i}(5)];end
Z4 = reorder_struct(Z4,listmap({'G(AGTC)C','C(AATC)G','G(ATTC)C','T(TTTC)A','T(CCTC)A','C(CCTC)G','T(GTTC)A'},Z4.context));
Z3 = reorder_struct(Z3,listmap({'G(GTCC)C','G(GTCT)C','T(TTCA)A','C(ATCT)G','C(CTCT)G','T(CTCA)A','C(CTCA)G'},Z3.context));
ZL3=reorder_struct(QL3,9);ZL3=rmfield(ZL3,{'min','max','label'});fs=fieldnames(ZL3);for fi=1:length(fs),f=fs{fi};ZL3.(f)=ZL3.(f)';end
ZL3.looplen=repmat(3,slength(ZL3),1);ZL3.looppos=repmat(3,slength(ZL3),1);ZL3.context=contextL3.name;
for i=1:slength(ZL3),ZL3.context{i}=[rc(ZL3.context{i}(4)) '(' ZL3.context{i}(1:3) ')' ZL3.context{i}(4)];end
ZL3all = reorder_struct(ZL3,listmap({'G(GTC)C';'G(ATC)C';'T(ATC)A';'T(GTC)A';'A(ATC)T';'G(CTC)C';'G(TTC)C';'A(TTC)T';'T(CTC)A';'A(CTC)T';'T(TTC)A';'A(GTC)T';...
                    'C(CTC)G';'C(GTC)G';'C(ATC)G';'C(TTC)G'},ZL3.context));
ZL3 = reorder_struct(ZL3,listmap({'G(GTC)C';'T(GTC)A';'A(GTC)T';'C(TTC)G'},ZL3.context));
save('FINAL.Qs2.v3a.1.mat','QL3','Q4','Q3','contextL3','context3','context4','Z3','Z4','ZL3','ZL3all');

% print results
pr(Z4);pr(Z3);pr(ZL3);pr(ZL3all)
%%%

% PLOT  Fig. 3f, S8

load('FINAL.Qs2.v3a.1.mat','Z3','Z4','ZL3','ZL3all');

figure(3),clf,fontsize=8;
subplot(2,2,1) % TETRALOOPS (Z3=pos3, Z4=pos4)
  dat=[Z4.relrate Z3.relrate]'; sd=[Z4.relsd Z3.relsd]'; context=[Z4.context;Z3.context];
  barweb(dat,dat-sd,dat+sd,0.8,{});colormap(jet);ff;set(gca,'fontsize',fontsize,'ytick',0:5:35);ylim([0 39]);
  xlabel('sequence of (loop) and closing basepair','fontsize',fontsize);  ylabel('relative mutation rate','fontsize',fontsize);
  line(xlim,[1 1],'linestyle','--','color',[0 0 0]); dd=dat';sdd=sd';
  for i=1:length(context),text(0.58+0.115*i+0.193*(i>7),dd(i)+sdd(i)+3.6,context{i},'fontsize',fontsize,'hor','cen','ver','bot','rot',72);end
subplot(2,2,2) % TRILOOPS (selected, for main figure)
  dat=[ZL3.relrate]'; sd=[ZL3.relsd]'; context=[ZL3.context];
  barweb(dat,dat-sd,dat+sd,0.8,{});colormap(jet);ff;set(gca,'fontsize',fontsize,'ytick',0:50:250);
  xlabel('sequence of (loop) and closing basepair','fontsize',fontsize);
  xlim([0.5 3.75]);ylim([0 260]);line(xlim,[1 1],'linestyle','--','color',[0 0 0]);
  for i=1:length(context),text(0.55+0.18*i,dat(i)+sd(i)+15,context{i},'fontsize',fontsize,'hor','cen','rot',72);end
subplot(2,2,3) % TRILOOPS (all--for supplementary info)
  dat=[ZL3all.relrate]'; sd=[ZL3all.relsd]'; context=[ZL3all.context];
  barweb(dat,dat-sd,dat+sd,0.8,{});colormap(jet);ff;set(gca,'fontsize',fontsize,'ytick',0:50:250);
  xlabel('sequence of (loop) and closing basepair','fontsize',fontsize); ylabel('relative mutation rate','fontsize',fontsize);
  xlim([0.58 1.4]);ylim([0 260]);  line(xlim,[1 1],'linestyle','--','color',[0 0 0]);
  for i=1:length(context),text(0.575+0.05*i,dat(i)+sd(i)+12,context{i},'fontsize',fontsize,'hor','cen','rot',70);end


% PLOT  Fig. S8

load('FINAL.Qs2.v3a.1.mat','QL3','contextL3');
ssbin=6; % this corresponds to stemstrength=13-14.  The actual stemstrength for this construct is 16-18... pretty close
context=contextL3.name; for i=1:length(context),context{i}=[rc(context{i}(4)) '(' context{i}(1:3) ')' context{i}(4)];end
Z=[]; Z.context=context; Z.relrate = QL3.relrate(ssbin,:)'; Z.relsd = QL3.relsd(ssbin,:)';
Z = reorder_struct(Z,grepm('\([GT]TC\)',Z.context)); Z = reorder_struct(Z,[3 7 1 5 4 8 2 6]); % reorder to match the figure

figure(9),clf
dat=Z.relrate; sd=Z.relsd; context=Z.context;
barweb(dat,dat-sd,dat+sd,0.8,{});colormap(jet);ff;set(gca,'fontsize',14,'linewidth',1.5); ylim([0 26.5]);
xlabel('sequence of (loop) and closing basepair','fontsize',14);  ylabel('relative mutation rate','fontsize',14);
xlim([0.55 1.45]); line(xlim,[1 1],'linestyle',':','color',[0 0 0]); dd=dat';sdd=sd';
for i=1:length(context),text(0.56+0.0995*i,dd(i)+sdd(i)+2.0,context{i},'fontsize',17,'hor','cen','ver','bot','rot',55);end
set(gcf,'papersize',[10 7.5],'paperposition',[0.2 0.2 10-0.4 7.5-0.4]);
print_to_file('finalfigs/sampledata/v1/fig.S8.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WXS hotspots
% from MC3

load('all_12086_WXS_hotspots.2.3.mat','S');
f=grep('relrate|qidx|column',fieldnames(S));S=rename_fields(S,f,prefix(f,'old_'));

% annotate with status in latest CGC driver list
C=mals('Census_allTue_Mar_5_12_59_32_2019.tsv');
C=reorder_struct(C,grepm('Mis|N|S',C.MutationTypes)); % --> 348 mutation-driven genes
S=rmfield(S,'driver'); S.cgc = ismember(S.gene,C.GeneSymbol);

save('all_12086_WXS_hotspots.2.4.mat','S');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% QUANTITATIVE MODELING

tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
X.site.ss = X.site.nbp+2*X.site.ngc; ssmin = [0 4:22]; X.site.ssbin=min(max(1,X.site.ss-3),length(ssmin));
X.mut=reorder_struct(X.mut,~isnan(X.mut.sidx));
X.mut.muttype = nan(slength(X.mut),1);
X.mut.muttype((X.mut.ref==2 & X.mut.alt==3)|(X.mut.ref==3 & X.mut.alt==2))=1; % C->G
X.mut.muttype((X.mut.ref==2 & X.mut.alt==1)|(X.mut.ref==3 & X.mut.alt==4))=2; % C->A
X.mut.muttype((X.mut.ref==2 & X.mut.alt==4)|(X.mut.ref==3 & X.mut.alt==1))=3; % C->T

% COHORT definition
C=[];C.name={}; C.puse=[]; pat=X.pat; pat.idx=(1:slength(pat))';
pat = reorder_struct(pat,pat.msupe_neg); % exclude patients with >10% of MSI, smoking, UV, POLE, or ESO.
pat=sort_struct(pat,'frac_apobec',-1);
C.name{end+1,1}='APOBEC>=10%'; C.puse(end+1,:)=ismember((1:slength(X.pat)),pat.idx(pat.frac_apobec>=0.1));
C.npat = sum(C.puse,2);
for ci=1:length(C.name),fprintf('%d/%d ',ci,length(C.name));
  C.muse{ci,1} = C.puse(ci,X.mut.pat_idx)==1; C.nmut(ci,1) = sum(C.muse{ci});
end,fprintf('\n');pr(C)

% TABULATE 
tic
maxlooplen=11; maxssbin=20; nq = round(1.1*sum(3:maxlooplen)*(4.^4)); % overestimte size of array to allocate
mo = keep_fields(X.mut,{'sidx','muttype'});
So = keep_fields(X.site,{'zone','looplen','looppos','minus2','minus1','plus1','plus2','ssbin'});
for ci=1:length(C.name),fprintf('cohort %d/%d:',ci,length(C.name));
  m = reorder_struct(mo,C.muse{ci}); S = So; S.ct = hist2d_fast_wrapper(m.sidx,m.muttype,1,slength(S),1,3);
  S = reorder_struct(rmfield(S,'zone'),S.zone<3); % exclude exons+spliceflank=20
  % TABULATE n and N across (looplen=3-maxlooplen, looppos=1-looplen, minus2=1-4, minus1=1-4, plus1=1-4, plus2=1-4, ssbin=1-maxssbin, muttype=1-3)
  N = nan(9,maxlooplen,4,4,4,4,maxssbin); n = nan(9,maxlooplen,4,4,4,4,maxssbin,3); fprintf(' looptype');
  for looplen=3:maxlooplen, S1 = reorder_struct(rmfield(S,'looplen'),S.looplen==looplen);
    for looppos=1:looplen, S2 = reorder_struct(rmfield(S1,'looppos'),S1.looppos==looppos); fprintf(' %d/%d',looppos,looplen);
      for minus2=1:4, S3=reorder_struct(rmfield(S2,'minus2'),S2.minus2==minus2);
        for minus1=1:4, S4=reorder_struct(rmfield(S3,'minus1'),S3.minus1==minus1);
          for plus1=1:4, S5=reorder_struct(rmfield(S4,'plus1'),S4.plus1==plus1);
            for plus2=1:4, S6=reorder_struct(rmfield(S5,'plus2'),S5.plus2==plus2);
              for ssbin=1:maxssbin, S7=reorder_struct(rmfield(S6,'ssbin'),S6.ssbin==ssbin);
                N(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin) = slength(S7);
                n(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin,:) = sum(S7.ct,1);
  end,end,end,end,end,end,end,fprintf('\n');
  % SECOND LEVEL OF TABULATION: make list Q of stemstrength series for each looptype+context
  Q=[];qi=0; fs={'looplen','looppos','minus2','minus1','plus1','plus2'}; for fi=1:length(fs),f=fs{fi};Q.(f)=nan(nq,1); end; Q.N=nan(nq,maxssbin); Q.n3=nan(nq,maxssbin,3);
  for looplen=3:maxlooplen, N1=squeeze(N(looplen-2,:,:,:,:,:,:)); n1=squeeze(n(looplen-2,:,:,:,:,:,:,:));
    for looppos=1:looplen, N2=squeeze(N1(looppos,:,:,:,:,:)); n2=squeeze(n1(looppos,:,:,:,:,:,:));
      for minus2=0:4, if minus2==0, N3=squeeze(sum(N2,1)); n3=squeeze(sum(n2,1)); else N3=squeeze(N2(minus2,:,:,:,:)); n3=squeeze(n2(minus2,:,:,:,:,:)); end
        for minus1=0:4, if minus1==0, N4=squeeze(sum(N3,1)); n4=squeeze(sum(n3,1)); else N4=squeeze(N3(minus1,:,:,:)); n4=squeeze(n3(minus1,:,:,:,:)); end
          for plus1=0:4, if plus1==0, N5=squeeze(sum(N4,1)); n5=squeeze(sum(n4,1)); else N5=squeeze(N4(plus1,:,:)); n5=squeeze(n4(plus1,:,:,:)); end
            for plus2=0:4, if plus2==0, N6=squeeze(sum(N5,1)); n6=squeeze(sum(n5,1)); else N6=squeeze(N5(plus2,:)); n6=squeeze(n5(plus2,:,:)); end
              qi=qi+1; for fi=1:length(fs),f=fs{fi};Q.(f)(qi)=eval(f);end; Q.N(qi,:)=N6; Q.n3(qi,:,:)=n6;
  end,end,end,end,end,end, Q=reorder_struct(Q,1:qi); Q.n=sum(Q.n3,3); [Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N); C.Q{ci,1}=Q;
end,toc % 10 min per cohort
% REFORMAT AND SAVE
Z=[]; Z.cohort=rmfield(C,'Q'); Z.stem.ssmin=as_column(ssmin); Z.muttype.name={'C->G';'C->A';'C->T'};
fs={'n3','n','rate','sd'}; Z.loop=rmfield(C.Q{1},fs); for fi=1:length(fs),f=fs{fi};for ci=1:slength(C); Z.loop.(f)(:,:,ci,:)=C.Q{ci}.(f);end,end
save('FINAL_cohorts_vs_looptypes.tabulated.1.0.mat','Z');
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CURVE FITTING FOR EACH SERIES
% Note, here "fo" and "sso" are referred to as "rinit" and "xa"

tmp=load('FINAL_cohorts_vs_looptypes.tabulated.1.0.mat','Z');X=tmp.Z; % ~5 sec to load
minmut=8; xas=[0:0.5:50 51:1:100 102:2:200]; mhps=[0.30:0.01:0.60];
ssmin = X.stem.ssmin; nss=length(ssmin); x = ssmin; measurepoint=find(x==20);
for ci=1:slength(X.cohort),fprintf('Cohort %d/%d\n',ci,length(X.cohort.name)); npat=X.cohort.npat(ci);
  Q=X.loop; fs={'n3','n','rate','sd'};for fi=1:length(fs),f=fs{fi};Q.(f)=Q.(f)(:,:,ci,:); end;
  W=[]; W.mhp=as_column(mhps); Q.rate=1e6*Q.rate/npat;Q.sd=1e6*Q.sd/npat;
  for wi=1:slength(W),fprintf('%d/%d ',wi,length(W.mhp)); mhp=W.mhp(wi); r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;
    rhp = nan(length(x),length(xas)); for xai=1:length(xas), rhp(:,xai) = r0*2.^(mhp*(x-xas(xai))); end
    for qi=slength(Q):-1:1
      robs = as_column(Q.rate(qi,:)); sdobs = as_column(Q.sd(qi,:)); robs(Q.n(qi,:)<minmut)=nan;
      rinit = robs(1); Q.rinit(qi,1)=rinit; rexp = rinit+rhp; yobs = log10(0.0001+robs); yexp = log10(0.0001+rexp);
      yerr = nansum(bsxfun(@minus,yexp,yobs).^2,1); [err xai] = min(yerr);
      Q.err(qi,1)=err; Q.xa(qi,1)=xas(xai); Q.rexp(qi,:) = rexp(:,xai);
      lastdata = find(~isnan(robs),1,'last'); if isempty(lastdata), lastdata=1; end % cap predictions where data runs out
      Q.rexp(qi,lastdata+1:end)=Q.rexp(qi,lastdata); Q.relrate_exp(qi,:) = Q.rexp(qi,:)/r0; Q.rr = Q.relrate_exp(:,measurepoint);
      Q.rexp(qi,lastdata+1:end)=nan; Q.lastdata(qi,1)=lastdata;
    end,W.err(wi,1)=sum(Q.err); W.Q{wi,1}=Q;
  end,fprintf('\n'); [tmp wi] = min(W.err); X.cohort.mhp(ci,1)=W.mhp(wi); X.cohort.Q{ci,1}=W.Q{wi};
end % ~10 sec per cohort
save('FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SHOW FITS
% on triangular arrangement of subplots

load('FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X');
ci=1; Q=X.cohort.Q{ci}; npat=X.cohort.npat(ci);r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;
figure(3),clf;nrow=9;ncol=11;minmut=10; x=X.stem.ssmin;
for looplen=3:11, for looppos=1:looplen, qi=find(Q.looplen==looplen & Q.looppos==looppos & Q.minus2==0 & Q.minus1==0 & Q.plus1==0 & Q.plus2==0);
  tpos = looppos-(looplen/2); row=looplen-2; col=(ncol/2)+tpos; xgap=0.018; ygap=0.028; width=(1/(0.7+ncol)); height=(1/(0.5+nrow));
  bottom = 1-row*height; left=0.8*xgap+(col-0.5)*(width); subplot('position',[left bottom width-xgap height-ygap]);
  robs = as_column(Q.rate(qi,:)); robs(Q.n(qi,:)<minmut)=nan; rinit = robs(1); sdobs = as_column(Q.sd(qi,:)); rexp = Q.rexp(qi,:);
  yobs=log10(robs);yinit=yobs(1);ymax=nanmax(yobs);yexp=log10(rexp);ylo=log10(max(0.01,robs-1.96*sdobs));yhi=log10(max(0.01,robs+1.96*sdobs));yopthresh=log10(r0*4);
  red=sigmoid(Q.rr(qi)/4,1,2.3); optclr=0.95*[1 1-red 1-red];
  cla,hold on, plot(x,yobs);plot(x,yexp);ff;xlim([0 max(x)+1]);ylim(log10([10 10000]));title([num2str(looppos) '/' num2str(looplen)],'fontsize',14);
  rectangle('position',[min(xlim) min(ylim) diff(xlim) abs(yopthresh-min(ylim))],'linestyle','none','facecolor',[0.8 0.8 0.8]);
  rectangle('position',[min(xlim) yopthresh diff(xlim) abs(max(ylim)-yopthresh)],'linestyle','none','facecolor',optclr);
  plot(x,yexp,'-','linewidth',2,'color',[0.4 0.4 1]); plot(x,yobs,'k.','markersize',12);
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  if looppos==1, ylabel('muts/Mb'); yt=[10 100 1000 10000]; set(gca,'ytick',log10(yt),'yticklabel',yt); else set(gca,'ytick',[]); end
  if row==nrow, xlabel('stem strength');set(gca,'xtick',[0 10 20],'xticklabel',[0 10 20]);else set(gca,'xtick',[]); end
end,end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SHOW GRIDS OF FITS FOR plus1/minus1

load('FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X');
ci=1; Q=X.cohort.Q{ci}; npat=X.cohort.npat(ci);r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;

looplen=3;looppos=3;
figure(4),clf;nrow=4;ncol=4;minmut=5;x=X.stem.ssmin;acgt='ACGT';
for minus1=1:4, for plus1=1:4, qi = find(Q.looplen==looplen & Q.looppos==looppos & Q.minus2==0 & Q.minus1==minus1 & Q.plus1==plus1 & Q.plus2==0);
  col=minus1;row=plus1;subplot(nrow,ncol,(row-1)*ncol+col); ttl=[acgt(5-plus1) '(' acgt(minus1) 'TC)' acgt(plus1)];
  robs = as_column(Q.rate(qi,:)); robs(Q.n(qi,:)<minmut)=nan; rinit = robs(1); sdobs = as_column(Q.sd(qi,:)); rexp = Q.rexp(qi,:);
  yobs=log10(robs);yinit=yobs(1);ymax=nanmax(yobs);yexp=log10(rexp);ylo=log10(max(0.01,robs-1.96*sdobs));yhi=log10(max(0.01,robs+1.96*sdobs));yopthresh=log10(r0*4);
  red=sigmoid(Q.rr(qi)/4,1,2.3); optclr=0.95*[1 1-red 1-red];
  cla,hold on, plot(x,yobs);plot(x,yexp);ff;xlim([0 max(x)+1]);ylim(log10([3 10000]));title(ttl,'fontsize',18);set(gca,'fontsize',14);
  rectangle('position',[min(xlim) min(ylim) diff(xlim) abs(yopthresh-min(ylim))],'linestyle','none','facecolor',[0.8 0.8 0.8]);
  rectangle('position',[min(xlim) yopthresh diff(xlim) abs(max(ylim)-yopthresh)],'linestyle','none','facecolor',optclr);
  plot(x,yexp,'-','linewidth',2,'color',[0.4 0.4 1]); plot(x,yobs,'k.','markersize',12);
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  if col==1, ylabel('muts/Mb','fontsize',16); yt=[10 100 1000 10000]; set(gca,'ytick',log10(yt),'yticklabel',yt); else set(gca,'ytick',[]); end
  if row==nrow, xlabel('stem strength','fontsize',16);set(gca,'xtick',[0 10 20],'xticklabel',[0 10 20]);else set(gca,'xtick',[]); end
end,end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% APPLY MODEL TO WXS HOTSPOTS TO DISTINGUISH
% PASSENGERS VS. DRIVERS

load('all_12086_WXS_hotspots.2.4.mat','S');

MODEL=load('FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X');MODEL=MODEL.X;
ci=1; Q=MODEL.cohort.Q{ci}; Q = reorder_struct(rmfield(Q,{'minus2','plus2'}),Q.minus2==0 & Q.plus2==0);

S.qidx1 = multimap(S,Q,{'looplen','looppos','minus1','plus1'});
S0=S;S0.minus1(:)=0;S0.plus1(:)=0; S.qidx0 = multimap(S0,Q,{'looplen','looppos','minus1','plus1'});
nssbin=20; S.ssbin=min(nssbin,max(1,S.dG3-3));
S.relrate_exp0 = nan(slength(S),1);
S.relrate_exp1 = nan(slength(S),1);
for i=1:slength(S),x=S.ssbin(i);
  S.relrate_exp0(i) = Q.relrate_exp(S.qidx0(i),x);
  S.relrate_exp1(i) = Q.relrate_exp(S.qidx1(i),x);
end
S.relrate_exp = S.relrate_exp1; idx=find(isnan(S.relrate_exp));S.relrate_exp(idx)=S.relrate_exp0(idx);

prn(S,{'minus1','plus1','dG3','looplen','looppos','old_relrate_exp','relrate_exp','driver','gene','flag'},1:100)

save('FINAL.all_12086_WXS_hotspots.3.0.mat','S');
%%%




%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIGURE 4

load('FINAL.all_12086_WXS_hotspots.3.0.mat','S');
n_show=100; S=reorder_struct(S,1:n_show);

figure(2);clf,hold on;randseed=3;randinit(randseed);
x = logjit(S.relrate_exp);y = logjit(S.ct);
clr = repmat([0.4 0.4 0.4],slength(S),1);
idx=find(S.cgc); clr(idx,:)=repmat([0 0 1],length(idx),1);
idx=find(~S.cgc & S.relrate_exp<4); clr(idx,:)=repmat([1 0 0],length(idx),1);
scatter(x,y,50,clr,'filled');ff;ylim([0.70 2.26]);xlim([0.12 2.2]);set(gca,'fontsize',20);
xlabel('substrate optimality','fontsize',30);ylabel('# mutated patients','fontsize',30);
xpts=[1 4 6 10 15 20 30 50 70 100 150]; set(gca,'xtick',log10(1.5+xpts),'xticklabel',prefix(num2cellstr(xpts),' '));
ypts=[  4 6 10 15 20 30 50 70 100 150]; set(gca,'ytick',log10(1.5+ypts),'yticklabel',ypts);
S.showname = ((S.ct>=7 & S.relrate_exp<4) | (S.ct>=6 & S.relrate_exp>=4)); fsz=repmat(10,slength(S),1);
tested = ismember(S.gene,{'TBC1D12','MB21D2','C3orf70','MROH2B','RARS2','NUP93'}); fsz(tested)=14;
idx=find(tested);scatter(x(idx),y(idx),100,clr(idx,:),'filled'); scatter(x(idx),y(idx),220,[0 0 0],'linewidth',1.5); S.showname(idx)=true;
feature = (~S.cgc & S.relrate_exp<4 & S.ct>=7); idx=find(feature);scatter(x(idx),y(idx),50,clr(idx,:),'filled');S.showname(idx)=true; fsz(idx)=14;
P=[];P.xpad=1; P.ypad=1;
idx=find(S.showname); textfit(x(idx)+(0.010+0.004*tested(idx))*diff(xlim),y(idx),S.gene(idx),'fontsize',fsz(idx),'fontcolor',clr(idx,:),P);
% --> NOTE: need to manually enlarge figure window and then repeat rendering, to get textit spacing right, *then* write to PDF.
set(gcf,'papersize',[12 10],'paperposition',[0.2 0.2 12-0.4 10-0.4]);
print_to_file('finalfigs/sampledata/v1/fig.4.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% APPLY MODEL GENOMEWIDE
%
% Annotate each TpC site withrelrate_exp
% --> in cases where there is enough data to use minus1,plus1, use those to refine rate.
% --> otherwise, just use looplen/looppos

tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec

MODEL=load('FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X'); MODEL=MODEL.X;
ci=1; Q=MODEL.cohort.Q{ci}; Q = reorder_struct(rmfield(Q,{'minus2','plus2'}),Q.minus2==0 & Q.plus2==0);

X.site.ss=X.site.nbp+2*X.site.ngc; nssbin=slength(MODEL.stem); X.site.ssbin=min(nssbin,max(1,X.site.ss-3));
z=nan(slength(X.site),1,'single'); X.site.relrate_exp0=z; X.site.relrate_exp1=z; tic
for looplen=3:11, i1=find(X.site.looplen==looplen);q1=find(Q.looplen==looplen);
  for looppos=1:looplen, i2=i1(X.site.looppos(i1)==looppos);q2=q1(Q.looppos(q1)==looppos); fprintf('%d/%d ',looppos,looplen);
    for ssbin=1:nssbin, i3=i2(X.site.ssbin(i2)==ssbin);q3=q2(Q.minus1(q2)==0 & Q.plus1(q2)==0); x=ssbin;
      y=q3; X.site.relrate_exp0(i3)=Q.relrate_exp(y,x);
      for minus1=1:4, i4=i3(X.site.minus1(i3)==minus1);q4=q2(Q.minus1(q2)==minus1);
        for plus1=1:4, i5=i4(X.site.plus1(i4)==plus1);q5=q4(Q.plus1(q4)==plus1);
          y=q5; X.site.relrate_exp1(i5)=Q.relrate_exp(y,x);
end,end,end,end,end, fprintf('\n'); toc % ~15 min
X.site.relrate_exp = X.site.relrate_exp1; idx=find(isnan(X.site.relrate_exp1)); X.site.relrate_exp(idx)=X.site.relrate_exp0(idx);
X.site=rmfield(X.site,{'ssbin','relrate_exp0','relrate_exp1'});
save('FINAL_DATASET.with_TpCs.with_model.v1.0.mat','X','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIGURE S4b (% mutations at optimal hairpins)

load('FINAL_DATASET.with_TpCs.with_model.v1.1.mat','X');
pat=X.pat;

figure(1),clf,hold on
  x = 100*pat.frac_apobec;
  y = 100*sum(pat.frrbin(:,5:end),2);
  clr = repmat([0.8 0.8 0.8],slength(pat),1); sz=10*ones(slength(pat),1);
  van = (pat.vanilla); clr(van,:)=repmat([0.2 0.7 0.2],sum(van),1);sz(van)=30;
  apo = (pat.apobec); clr(apo,:)=repmat([0 0 0],sum(apo),1);sz(apo)=30;
  a3b = (pat.apobec_most_a3b);clr(a3b,:)=repmat([0 0 1],sum(a3b),1); sz(a3b)=50;
  a3a = (pat.apobec_most_a3a);clr(a3a,:)=repmat([1 0 0],sum(a3a),1); sz(a3a)=50;
  abd = (apo & pat.ABdel==1);
  scatter(x,y,sz,clr,'filled');
  scatter(x(van),y(van),sz(van),clr(van,:),'filled');
  scatter(x(apo),y(apo),sz(apo),clr(apo,:),'filled');
  scatter(x(a3b),y(a3b),sz(a3b),clr(a3b,:),'filled');
  scatter(x(a3a),y(a3a),sz(a3a),clr(a3a,:),'filled');
  scatter(x(abd),y(abd),130,[1 0.7 0],'filled');scatter(x(abd),y(abd),150,[0 0 0]);scatter(x(abd),y(abd),sz(abd),clr(abd,:),'filled');
  line([10 10],ylim,'linestyle','--','color',[0 0 0]);ff;set(gca,'fontsize',20);ylim([0 1.8]);
  xlabel('% APOBEC mutations','fontsize',30);ylabel('% at optimal hairpins','fontsize',30);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% VALIDATION OF MODELING

% (1) Predict WGS coding rates from WGS noncoding model
% (2) Predict WXS coding rates from WGS noncoding model

load('FINAL_DATASET.with_TpCs.with_model.v1.0.mat','X');
puse=(X.pat.apobec); muse=puse(X.mut.pat_idx); % msupe- apo+
Sx = X.site; Sx.ct_apobec_wgs = histc(X.mut.sidx(muse),1:slength(Sx));
Sx.idx_orig = (1:slength(Sx))'; Sx=off(Sx,'idx_orig'); Sx = reorder_struct(rmfield(Sx,'zone'),Sx.zone==3);
save('FINAL_DATASET.exonic_TpCs.with_model.v1.0.mat','Sx');
%%%
load('FINAL_DATASET.exonic_TpCs.with_model.v1.0.mat','Sx');
load('MC3.exome_TpC_sites.spliceflank10.1.4.Gflipped.compact.1.0.mat','M');
M.pat.nmf_norm = bsxfun(@rdivide,M.pat.nmf,nansum(M.pat.nmf,2));
M.pat.frac_apobec = M.pat.nmf_norm(:,9);
M.pat.frac_eso = M.pat.nmf_norm(:,1);
M.pat.frac_smoking = M.pat.nmf_norm(:,2);
M.pat.frac_pole = M.pat.nmf_norm(:,3);
M.pat.frac_uv = M.pat.nmf_norm(:,5);
M.pat.frac_msi = sum(M.pat.nmf_norm(:,[7 11 12]),2);
M.pat.msupe_neg = (M.pat.frac_msi<0.1 & M.pat.frac_smoking<0.1 & M.pat.frac_uv<0.1 & M.pat.frac_pole<0.1 & M.pat.frac_eso<0.1);
M.pat.apobec = (M.pat.frac_apobec>=0.1 & M.pat.msupe_neg);
M.mut.tpc = ismember(M.mut.context65,[29:33 37 41 45]);
m = reorder_struct(M.mut,M.mut.tpc & M.pat.apobec(M.mut.pat_idx));
tic; m.Sxidx = multimap(m,Sx,{'chr','pos'}); toc % ~10 sec
Sx.ct_apobec_wxs = histc(m.Sxidx,1:slength(Sx)); % <1 sec

save('validation.v1.0.mat','Sx');
%%%

load('validation.v1.0.mat','Sx');

V={};

bin=[]; bin.min=[0.2:0.1:2 3:1:10 12 15 20 30 50]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (Sx.relrate_exp>=bin.min(i) & Sx.relrate_exp<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(Sx.ct_apobec_wgs(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(Sx.ct_apobec_wgs); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;

V{1}=bin;

bin=[]; bin.min=[0.2:0.1:3 3.5:0.5:10 12:1:20 25 30 40 100 150]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (Sx.relrate_exp>=bin.min(i) & Sx.relrate_exp<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(Sx.ct_apobec_wxs(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(Sx.ct_apobec_wxs); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;

V{2}=bin;

save('validation.tabulated.v1.0.mat','V');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next:
% (3) split patients in half, predict one half from the other

% REPEAT TABULATION AND CURVE FITTING FOR THESE

tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
X.site.ss = X.site.nbp+2*X.site.ngc; ssmin = [0 4:22]; X.site.ssbin=min(max(1,X.site.ss-3),length(ssmin));
X.mut=reorder_struct(X.mut,~isnan(X.mut.sidx));
X.mut.muttype = nan(slength(X.mut),1);
X.mut.muttype((X.mut.ref==2 & X.mut.alt==3)|(X.mut.ref==3 & X.mut.alt==2))=1; % C->G
X.mut.muttype((X.mut.ref==2 & X.mut.alt==1)|(X.mut.ref==3 & X.mut.alt==4))=2; % C->A
X.mut.muttype((X.mut.ref==2 & X.mut.alt==4)|(X.mut.ref==3 & X.mut.alt==1))=3; % C->T

% COHORT definition

pat=X.pat; pat.idx=(1:slength(pat))';
pat = reorder_struct(pat,pat.apobec); % msupe- (<10%) apo+ (>=10%)
pat = sort_struct(pat,'nmut',-1);
pat.odd = mod((1:slength(pat))',2)==1;
pat.even = mod((1:slength(pat))',2)==0;

C=[];C.name={};
ci=1; C.name{ci,1}='APOBEC>=10% (odd patients)'; C.puse(ci,:)=ismember((1:slength(X.pat)),pat.idx(pat.odd)); C.muse(ci,:)=C.puse(ci,X.mut.pat_idx);
C.npat = sum(C.puse,2); C.nmut = sum(C.muse,2); pr(C)

% TABULATE
tic
maxlooplen=11; maxssbin=20; nq = round(1.1*sum(3:maxlooplen)*(4.^4));
mo = keep_fields(X.mut,{'sidx','muttype'});
So = keep_fields(X.site,{'zone','looplen','looppos','minus2','minus1','plus1','plus2','ssbin'});
for ci=1:length(C.name),fprintf('cohort %d/%d:',ci,length(C.name));
  m = reorder_struct(mo,C.muse(ci,:)); S = So; S.ct = hist2d_fast_wrapper(m.sidx,m.muttype,1,slength(S),1,3);
  S = reorder_struct(rmfield(S,'zone'),S.zone<3); % exclude exons+spliceflank=20
  % TABULATE n and N across (looplen=3-maxlooplen, looppos=1-looplen, minus2=1-4, minus1=1-4, plus1=1-4, plus2=1-4, ssbin=1-maxssbin, muttype=1-3)
  N = nan(9,maxlooplen,4,4,4,4,maxssbin); n = nan(9,maxlooplen,4,4,4,4,maxssbin,3); fprintf(' looptype');
  for looplen=3:maxlooplen, S1 = reorder_struct(rmfield(S,'looplen'),S.looplen==looplen);
    for looppos=1:looplen, S2 = reorder_struct(rmfield(S1,'looppos'),S1.looppos==looppos); fprintf(' %d/%d',looppos,looplen);
      for minus2=1:4, S3=reorder_struct(rmfield(S2,'minus2'),S2.minus2==minus2);
        for minus1=1:4, S4=reorder_struct(rmfield(S3,'minus1'),S3.minus1==minus1);
          for plus1=1:4, S5=reorder_struct(rmfield(S4,'plus1'),S4.plus1==plus1);
            for plus2=1:4, S6=reorder_struct(rmfield(S5,'plus2'),S5.plus2==plus2);
              for ssbin=1:maxssbin, S7=reorder_struct(rmfield(S6,'ssbin'),S6.ssbin==ssbin);
                N(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin) = slength(S7);
                n(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin,:) = sum(S7.ct,1);
  end,end,end,end,end,end,end,fprintf('\n');
  % SECOND LEVEL OF TABULATION: make list Q of stemstrength series for each looptype+context
  Q=[];qi=0; fs={'looplen','looppos','minus2','minus1','plus1','plus2'}; for fi=1:length(fs),f=fs{fi};Q.(f)=nan(nq,1); end; Q.N=nan(nq,maxssbin); Q.n3=nan(nq,maxssbin,3);
  for looplen=3:maxlooplen, N1=squeeze(N(looplen-2,:,:,:,:,:,:)); n1=squeeze(n(looplen-2,:,:,:,:,:,:,:));
    for looppos=1:looplen, N2=squeeze(N1(looppos,:,:,:,:,:)); n2=squeeze(n1(looppos,:,:,:,:,:,:));
      for minus2=0:4, if minus2==0, N3=squeeze(sum(N2,1)); n3=squeeze(sum(n2,1)); else N3=squeeze(N2(minus2,:,:,:,:)); n3=squeeze(n2(minus2,:,:,:,:,:)); end
        for minus1=0:4, if minus1==0, N4=squeeze(sum(N3,1)); n4=squeeze(sum(n3,1)); else N4=squeeze(N3(minus1,:,:,:)); n4=squeeze(n3(minus1,:,:,:,:)); end
          for plus1=0:4, if plus1==0, N5=squeeze(sum(N4,1)); n5=squeeze(sum(n4,1)); else N5=squeeze(N4(plus1,:,:)); n5=squeeze(n4(plus1,:,:,:)); end
            for plus2=0:4, if plus2==0, N6=squeeze(sum(N5,1)); n6=squeeze(sum(n5,1)); else N6=squeeze(N5(plus2,:)); n6=squeeze(n5(plus2,:,:)); end
              qi=qi+1; for fi=1:length(fs),f=fs{fi};Q.(f)(qi)=eval(f);end; Q.N(qi,:)=N6; Q.n3(qi,:,:)=n6;
  end,end,end,end,end,end, Q=reorder_struct(Q,1:qi); Q.n=sum(Q.n3,3); [Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N); C.Q{ci,1}=Q;
end,toc % 10 min per cohort
% REFORMAT AND SAVE
Z=[]; Z.cohort=rmfield(C,'Q'); Z.stem.ssmin=as_column(ssmin); Z.muttype.name={'C->G';'C->A';'C->T'};
fs={'n3','n','rate','sd'}; Z.loop=rmfield(C.Q{1},fs); for fi=1:length(fs),f=fs{fi};for ci=1:slength(C); Z.loop.(f)(:,:,ci,:)=C.Q{ci}.(f);end,end
save('validation_cohorts_vs_looptypes.tabulated.2.0.mat','Z');

% CURVE FITTING
tmp=load('validation_cohorts_vs_looptypes.tabulated.2.0.mat','Z');X=tmp.Z; % ~5 sec to load
minmut=8; xas=[0:0.5:50 51:1:100 102:2:200]; mhps=[0.30:0.01:0.60]; 
ssmin = X.stem.ssmin; nss=length(ssmin); x = ssmin; measurepoint=find(x==20);
for ci=1:slength(X.cohort),fprintf('Cohort %d/%d\n',ci,length(X.cohort.name)); npat=X.cohort.npat(ci);
  Q=X.loop; fs={'n3','n','rate','sd'};for fi=1:length(fs),f=fs{fi};Q.(f)=Q.(f)(:,:,ci,:); end;
  W=[]; W.mhp=as_column(mhps); Q.rate=1e6*Q.rate/npat;Q.sd=1e6*Q.sd/npat;
  for wi=1:slength(W),fprintf('%d/%d ',wi,length(W.mhp)); mhp=W.mhp(wi); r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;
    rhp = nan(length(x),length(xas)); for xai=1:length(xas), rhp(:,xai) = r0*2.^(mhp*(x-xas(xai))); end
    for qi=slength(Q):-1:1
      robs = as_column(Q.rate(qi,:)); sdobs = as_column(Q.sd(qi,:)); robs(Q.n(qi,:)<minmut)=nan;
      rinit = robs(1); Q.rinit(qi,1)=rinit; rexp = rinit+rhp; yobs = log10(0.0001+robs); yexp = log10(0.0001+rexp);
      yerr = nansum(bsxfun(@minus,yexp,yobs).^2,1); [err xai] = min(yerr);
      Q.err(qi,1)=err; Q.xa(qi,1)=xas(xai); Q.rexp(qi,:) = rexp(:,xai);
      lastdata = find(~isnan(robs),1,'last'); if isempty(lastdata), lastdata=1; end % cap predictions where data runs out
      Q.rexp(qi,lastdata+1:end)=Q.rexp(qi,lastdata); Q.relrate_exp(qi,:) = Q.rexp(qi,:)/r0; Q.rr = Q.relrate_exp(:,measurepoint);
      Q.rexp(qi,lastdata+1:end)=nan; Q.lastdata(qi,1)=lastdata;
    end,W.err(wi,1)=sum(Q.err); W.Q{wi,1}=Q;
  end,fprintf('\n'); [tmp wi] = min(W.err); X.cohort.mhp(ci,1)=W.mhp(wi); X.cohort.Q{ci,1}=W.Q{wi};
end
save('validation_cohorts_vs_looptypes.tabulated.model.2.0.mat','X');

pr(X.cohort)
%%%

% Apply models genomewide
tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
MODEL=load('validation_cohorts_vs_looptypes.tabulated.model.2.0.mat','X'); MODEL=MODEL.X;
X.site.ss=X.site.nbp+2*X.site.ngc; nssbin=slength(MODEL.stem); X.site.ssbin=min(nssbin,max(1,X.site.ss-3));
X.site.relrate_exp_cohorts = nan(slength(X.site),slength(MODEL.cohort),'single');
for ci=1:slength(MODEL.cohort), fprintf('COHORT %d ',ci); tic
  Q=MODEL.cohort.Q{ci}; Q = reorder_struct(rmfield(Q,{'minus2','plus2'}),Q.minus2==0 & Q.plus2==0);
  rr0=nan(slength(X.site),1,'single'); rr1=rr0;
  for looplen=3:11, i1=find(X.site.looplen==looplen);q1=find(Q.looplen==looplen);
    for looppos=1:looplen, i2=i1(X.site.looppos(i1)==looppos);q2=q1(Q.looppos(q1)==looppos); fprintf('%d/%d ',looppos,looplen);
      for ssbin=1:nssbin, i3=i2(X.site.ssbin(i2)==ssbin);q3=q2(Q.minus1(q2)==0 & Q.plus1(q2)==0); x=ssbin;
        y=q3; rr0(i3)=Q.relrate_exp(y,x);
        for minus1=1:4, i4=i3(X.site.minus1(i3)==minus1);q4=q2(Q.minus1(q2)==minus1);
          for plus1=1:4, i5=i4(X.site.plus1(i4)==plus1);q5=q4(Q.plus1(q4)==plus1);
            y=q5; rr1(i5)=Q.relrate_exp(y,x);
  end,end,end,end,end, fprintf('\n'); toc % ~15 min
  X.site.relrate_exp_cohorts(:,ci) = rr1; idx=find(isnan(rr1)); X.site.relrate_exp_cohorts(idx,ci)=rr0(idx);
end % ~6 min per cohort
save('FINAL_DATASET.with_TpCs.with_validation_models.v2.0.mat','X','-v7.3');
%%%

% how well does model work? (predict even patients from odd patients)

load('FINAL_DATASET.with_TpCs.with_validation_models.v2.0.mat','X');
pat=X.pat; pat.idx=(1:slength(pat))';
pat = reorder_struct(pat,pat.apobec); % msupe- (<10%) apo+ (>=10%)
pat = sort_struct(pat,'nmut',-1);
pat.odd = mod((1:slength(pat))',2)==1;
pat.even = mod((1:slength(pat))',2)==0;
X.mut.odd = mod(X.mut.pos,2)==1;
X.mut.even = mod(X.mut.pos,2)==0;

X.site.relrate_exp_odd = X.site.relrate_exp_cohorts(:,1);  % rate predicted from odd patients in msupe- apo+ cohort
puse=ismember((1:slength(X.pat))',pat.idx(pat.even)); muse=puse(X.mut.pat_idx); % mutations in even patients
X.site.ct_apobec_wgs_even = histc(X.mut.sidx(muse),1:slength(X.site));

bin=[]; bin.min=[0.2:0.1:2 3:1:10 12 15 20 30 50 70 100]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (X.site.relrate_exp_odd>=bin.min(i) & X.site.relrate_exp_odd<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(X.site.ct_apobec_wgs_even(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(X.site.ct_apobec_wgs_even); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;
pr(bin)

V={}; V{3}=bin;
save('validation.tabulated.v2.0.mat','V');
%%%

load('validation.tabulated.v2.0.mat','V');

figure(1)
bin=V{3};
x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff
ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
xlabel('relative mutation rate (predicted, WGS noncoding, odd patients)','fontsize',20);
ylabel('relative mutation rate (observed, WGS noncoding, even patients)','fontsize',20);
line(xlim,xlim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (4) predict even sites from odd sites

% REPEAT TABULATION AND CURVE FITTING

tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
X.site.ss = X.site.nbp+2*X.site.ngc; ssmin = [0 4:22]; X.site.ssbin=min(max(1,X.site.ss-3),length(ssmin));
X.mut=reorder_struct(X.mut,~isnan(X.mut.sidx));
X.mut.muttype = nan(slength(X.mut),1);
X.mut.muttype((X.mut.ref==2 & X.mut.alt==3)|(X.mut.ref==3 & X.mut.alt==2))=1; % C->G
X.mut.muttype((X.mut.ref==2 & X.mut.alt==1)|(X.mut.ref==3 & X.mut.alt==4))=2; % C->A
X.mut.muttype((X.mut.ref==2 & X.mut.alt==4)|(X.mut.ref==3 & X.mut.alt==1))=3; % C->T

% COHORT definition

pat=X.pat; pat.idx=(1:slength(pat))';
pat = reorder_struct(pat,pat.apobec); % msupe- (<10%) apo+ (>=10%)
pat = sort_struct(pat,'nmut',-1);
X.site.odd = mod(X.site.pos,2)==1;
X.site.even = mod(X.site.pos,2)==0;

C=[];C.name={};
ci=1; C.name{ci,1}='APOBEC>=10%'; C.puse(ci,:)=ismember((1:slength(X.pat)),pat.idx); C.muse(ci,:)=C.puse(ci,X.mut.pat_idx);
C.npat = sum(C.puse,2); C.nmut = sum(C.muse,2); pr(C)

% TABULATE
tic
maxlooplen=11; maxssbin=20; nq = round(1.1*sum(3:maxlooplen)*(4.^4));
mo = keep_fields(X.mut,{'sidx','muttype'});
So = keep_fields(X.site,{'odd','zone','looplen','looppos','minus2','minus1','plus1','plus2','ssbin'});
for ci=1:length(C.name),fprintf('cohort %d/%d:',ci,length(C.name));
  m = reorder_struct(mo,C.muse(ci,:)); S = So; S.ct = hist2d_fast_wrapper(m.sidx,m.muttype,1,slength(S),1,3);
  S = reorder_struct(rmfield(S,'zone'),S.zone<3); % exclude exons+spliceflank=20
  S = reorder_struct(rmfield(S,'odd'),S.odd);  % predict on basis of odd sites only
  % TABULATE n and N across (looplen=3-maxlooplen, looppos=1-looplen, minus2=1-4, minus1=1-4, plus1=1-4, plus2=1-4, ssbin=1-maxssbin, muttype=1-3)
  N = nan(9,maxlooplen,4,4,4,4,maxssbin); n = nan(9,maxlooplen,4,4,4,4,maxssbin,3); fprintf(' looptype');
  for looplen=3:maxlooplen, S1 = reorder_struct(rmfield(S,'looplen'),S.looplen==looplen);
    for looppos=1:looplen, S2 = reorder_struct(rmfield(S1,'looppos'),S1.looppos==looppos); fprintf(' %d/%d',looppos,looplen);
      for minus2=1:4, S3=reorder_struct(rmfield(S2,'minus2'),S2.minus2==minus2);
        for minus1=1:4, S4=reorder_struct(rmfield(S3,'minus1'),S3.minus1==minus1);
          for plus1=1:4, S5=reorder_struct(rmfield(S4,'plus1'),S4.plus1==plus1);
            for plus2=1:4, S6=reorder_struct(rmfield(S5,'plus2'),S5.plus2==plus2);
              for ssbin=1:maxssbin, S7=reorder_struct(rmfield(S6,'ssbin'),S6.ssbin==ssbin);
                N(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin) = slength(S7);
                n(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin,:) = sum(S7.ct,1);
  end,end,end,end,end,end,end,fprintf('\n');
  % SECOND LEVEL OF TABULATION: make list Q of stemstrength series for each looptype+context
  Q=[];qi=0; fs={'looplen','looppos','minus2','minus1','plus1','plus2'}; for fi=1:length(fs),f=fs{fi};Q.(f)=nan(nq,1); end; Q.N=nan(nq,maxssbin); Q.n3=nan(nq,maxssbin,3);
  for looplen=3:maxlooplen, N1=squeeze(N(looplen-2,:,:,:,:,:,:)); n1=squeeze(n(looplen-2,:,:,:,:,:,:,:));
    for looppos=1:looplen, N2=squeeze(N1(looppos,:,:,:,:,:)); n2=squeeze(n1(looppos,:,:,:,:,:,:));
      for minus2=0:4, if minus2==0, N3=squeeze(sum(N2,1)); n3=squeeze(sum(n2,1)); else N3=squeeze(N2(minus2,:,:,:,:)); n3=squeeze(n2(minus2,:,:,:,:,:)); end
        for minus1=0:4, if minus1==0, N4=squeeze(sum(N3,1)); n4=squeeze(sum(n3,1)); else N4=squeeze(N3(minus1,:,:,:)); n4=squeeze(n3(minus1,:,:,:,:)); end
          for plus1=0:4, if plus1==0, N5=squeeze(sum(N4,1)); n5=squeeze(sum(n4,1)); else N5=squeeze(N4(plus1,:,:)); n5=squeeze(n4(plus1,:,:,:)); end
            for plus2=0:4, if plus2==0, N6=squeeze(sum(N5,1)); n6=squeeze(sum(n5,1)); else N6=squeeze(N5(plus2,:)); n6=squeeze(n5(plus2,:,:)); end
              qi=qi+1; for fi=1:length(fs),f=fs{fi};Q.(f)(qi)=eval(f);end; Q.N(qi,:)=N6; Q.n3(qi,:,:)=n6;
  end,end,end,end,end,end, Q=reorder_struct(Q,1:qi); Q.n=sum(Q.n3,3); [Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N); C.Q{ci,1}=Q;
end,toc % 10 min per cohort
% REFORMAT AND SAVE
Z=[]; Z.cohort=rmfield(C,'Q'); Z.stem.ssmin=as_column(ssmin); Z.muttype.name={'C->G';'C->A';'C->T'};
fs={'n3','n','rate','sd'}; Z.loop=rmfield(C.Q{1},fs); for fi=1:length(fs),f=fs{fi};for ci=1:slength(C); Z.loop.(f)(:,:,ci,:)=C.Q{ci}.(f);end,end
save('validation_cohorts_vs_looptypes.tabulated.3.0.mat','Z');

% CURVE FITTING
tmp=load('validation_cohorts_vs_looptypes.tabulated.3.0.mat','Z');X=tmp.Z; % ~5 sec to load
minmut=8; xas=[0:0.5:50 51:1:100 102:2:200]; mhps=[0.30:0.01:0.60];
ssmin = X.stem.ssmin; nss=length(ssmin); x = ssmin; measurepoint=find(x==20);
for ci=1:slength(X.cohort),fprintf('Cohort %d/%d\n',ci,length(X.cohort.name)); npat=X.cohort.npat(ci);
  Q=X.loop; fs={'n3','n','rate','sd'};for fi=1:length(fs),f=fs{fi};Q.(f)=Q.(f)(:,:,ci,:); end;
  W=[]; W.mhp=as_column(mhps); Q.rate=1e6*Q.rate/npat;Q.sd=1e6*Q.sd/npat;
  for wi=1:slength(W),fprintf('%d/%d ',wi,length(W.mhp)); mhp=W.mhp(wi); r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;
    rhp = nan(length(x),length(xas)); for xai=1:length(xas), rhp(:,xai) = r0*2.^(mhp*(x-xas(xai))); end
    for qi=slength(Q):-1:1
      robs = as_column(Q.rate(qi,:)); sdobs = as_column(Q.sd(qi,:)); robs(Q.n(qi,:)<minmut)=nan;
      rinit = robs(1); Q.rinit(qi,1)=rinit; rexp = rinit+rhp; yobs = log10(0.0001+robs); yexp = log10(0.0001+rexp);
      yerr = nansum(bsxfun(@minus,yexp,yobs).^2,1); [err xai] = min(yerr);
      Q.err(qi,1)=err; Q.xa(qi,1)=xas(xai); Q.rexp(qi,:) = rexp(:,xai);
      lastdata = find(~isnan(robs),1,'last'); if isempty(lastdata), lastdata=1; end % cap predictions where data runs out
      Q.rexp(qi,lastdata+1:end)=Q.rexp(qi,lastdata); Q.relrate_exp(qi,:) = Q.rexp(qi,:)/r0; Q.rr = Q.relrate_exp(:,measurepoint);
      Q.rexp(qi,lastdata+1:end)=nan; Q.lastdata(qi,1)=lastdata;
    end,W.err(wi,1)=sum(Q.err); W.Q{wi,1}=Q;
  end,fprintf('\n'); [tmp wi] = min(W.err); X.cohort.mhp(ci,1)=W.mhp(wi); X.cohort.Q{ci,1}=W.Q{wi};
end
save('validation_cohorts_vs_looptypes.tabulated.model.3.0.mat','X');

pr(X.cohort)

% Apply model genomewide
tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
MODEL=load('validation_cohorts_vs_looptypes.tabulated.model.3.0.mat','X'); MODEL=MODEL.X;
X.site.ss=X.site.nbp+2*X.site.ngc; nssbin=slength(MODEL.stem); X.site.ssbin=min(nssbin,max(1,X.site.ss-3));
X.site.relrate_exp_cohorts = nan(slength(X.site),slength(MODEL.cohort),'single');
for ci=1:slength(MODEL.cohort), fprintf('COHORT %d ',ci); tic
  Q=MODEL.cohort.Q{ci}; Q = reorder_struct(rmfield(Q,{'minus2','plus2'}),Q.minus2==0 & Q.plus2==0);
  rr0=nan(slength(X.site),1,'single'); rr1=rr0;
  for looplen=3:11, i1=find(X.site.looplen==looplen);q1=find(Q.looplen==looplen);
    for looppos=1:looplen, i2=i1(X.site.looppos(i1)==looppos);q2=q1(Q.looppos(q1)==looppos); fprintf('%d/%d ',looppos,looplen);
      for ssbin=1:nssbin, i3=i2(X.site.ssbin(i2)==ssbin);q3=q2(Q.minus1(q2)==0 & Q.plus1(q2)==0); x=ssbin;
        y=q3; rr0(i3)=Q.relrate_exp(y,x);
        for minus1=1:4, i4=i3(X.site.minus1(i3)==minus1);q4=q2(Q.minus1(q2)==minus1);
          for plus1=1:4, i5=i4(X.site.plus1(i4)==plus1);q5=q4(Q.plus1(q4)==plus1);
            y=q5; rr1(i5)=Q.relrate_exp(y,x);
  end,end,end,end,end, fprintf('\n'); toc % ~15 min
  X.site.relrate_exp_cohorts(:,ci) = rr1; idx=find(isnan(rr1)); X.site.relrate_exp_cohorts(idx,ci)=rr0(idx);
end % ~6 min per cohort
save('FINAL_DATASET.with_TpCs.with_validation_models.v3.0.mat','X','-v7.3');
%%%

% how well does model work? (predict even sites from odd sites)

load('FINAL_DATASET.with_TpCs.with_validation_models.v3.0.mat','X');
pat=X.pat; pat.idx=(1:slength(pat))';
pat = reorder_struct(pat,pat.apobec); % msupe- (<10%) apo+ (>=10%)
pat = sort_struct(pat,'nmut',-1);
X.site.odd = mod(X.site.pos,2)==1;
X.site.even = mod(X.site.pos,2)==0;

X.site.relrate_exp_odd = X.site.relrate_exp_cohorts(:,1);  % rate predicted from odd sites in msupe- apo+ cohort
puse=ismember((1:slength(X.pat))',pat.idx); muse=puse(X.mut.pat_idx);
X.site.ct_apobec_wgs = histc(X.mut.sidx(muse),1:slength(X.site));  % rate observed at all sites
S_even = reorder_struct(X.site,X.site.even);

bin=[]; bin.min=[0.2:0.1:2 3:1:10 12 15 20 30 50 70 100]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (S_even.relrate_exp_odd>=bin.min(i) & S_even.relrate_exp_odd<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(S_even.ct_apobec_wgs(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(S_even.ct_apobec_wgs); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;
pr(bin)

V={}; V{4}=bin;
save('validation.tabulated.v3.0.mat','V');
%%%

load('validation.tabulated.v3.0.mat','V');

figure(1)
bin=V{4};
x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff
ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
xlabel('relative mutation rate (predicted, WGS noncoding, odd sites)','fontsize',20);
ylabel('relative mutation rate (observed, WGS noncoding, even sites)','fontsize',20);
line(xlim,xlim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (5) make model based on a single patient's mutations

% name_________________  ttype nmut  frac_apobec
% Nik-Zainal560:PD13604a BRCA  93102 0.90
% Alexandrov507:PD4120a  BRCA  67364 0.94

tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
X.site.ss = X.site.nbp+2*X.site.ngc; ssmin = [0 4:22]; X.site.ssbin=min(max(1,X.site.ss-3),length(ssmin));
X.mut=reorder_struct(X.mut,~isnan(X.mut.sidx));
X.mut.muttype = nan(slength(X.mut),1);
X.mut.muttype((X.mut.ref==2 & X.mut.alt==3)|(X.mut.ref==3 & X.mut.alt==2))=1; % C->G
X.mut.muttype((X.mut.ref==2 & X.mut.alt==1)|(X.mut.ref==3 & X.mut.alt==4))=2; % C->A
X.mut.muttype((X.mut.ref==2 & X.mut.alt==4)|(X.mut.ref==3 & X.mut.alt==1))=3; % C->T

% COHORT definition

pat=X.pat; pat.idx=(1:slength(pat))';
pat = reorder_struct(pat,pat.apobec); % msupe- (<10%) apo+ (>=10%)
pat = sort_struct(pat,'nmut',-1);

C=[];C.name={};
ci=1; C.name{ci,1}='top APOBEC patient'; C.puse(ci,:)=ismember((1:slength(X.pat)),pat.idx(1)); C.muse(ci,:)=C.puse(ci,X.mut.pat_idx);
ci=2; C.name{ci,1}='second APOBEC patient'; C.puse(ci,:)=ismember((1:slength(X.pat)),pat.idx(2)); C.muse(ci,:)=C.puse(ci,X.mut.pat_idx);
C.npat = sum(C.puse,2); C.nmut = sum(C.muse,2); pr(C)

name_________________ puse  muse  npat nmut
top APOBEC patient    {...} {...} 1    81747
second APOBEC patient {...} {...} 1    61352

% TABULATE
tic
maxlooplen=11; maxssbin=20; nq = round(1.1*sum(3:maxlooplen)*(4.^4));
mo = keep_fields(X.mut,{'sidx','muttype'});
So = keep_fields(X.site,{'zone','looplen','looppos','minus2','minus1','plus1','plus2','ssbin'});
for ci=1:length(C.name),fprintf('cohort %d/%d:',ci,length(C.name));
  m = reorder_struct(mo,C.muse(ci,:)); S = So; S.ct = hist2d_fast_wrapper(m.sidx,m.muttype,1,slength(S),1,3);
  S = reorder_struct(rmfield(S,'zone'),S.zone<3); % exclude exons+spliceflank=20
  % TABULATE n and N across (looplen=3-maxlooplen, looppos=1-looplen, minus2=1-4, minus1=1-4, plus1=1-4, plus2=1-4, ssbin=1-maxssbin, muttype=1-3)
  N = nan(9,maxlooplen,4,4,4,4,maxssbin); n = nan(9,maxlooplen,4,4,4,4,maxssbin,3); fprintf(' looptype');
  for looplen=3:maxlooplen, S1 = reorder_struct(rmfield(S,'looplen'),S.looplen==looplen);
    for looppos=1:looplen, S2 = reorder_struct(rmfield(S1,'looppos'),S1.looppos==looppos); fprintf(' %d/%d',looppos,looplen);
      for minus2=1:4, S3=reorder_struct(rmfield(S2,'minus2'),S2.minus2==minus2);
        for minus1=1:4, S4=reorder_struct(rmfield(S3,'minus1'),S3.minus1==minus1);
          for plus1=1:4, S5=reorder_struct(rmfield(S4,'plus1'),S4.plus1==plus1);
            for plus2=1:4, S6=reorder_struct(rmfield(S5,'plus2'),S5.plus2==plus2);
              for ssbin=1:maxssbin, S7=reorder_struct(rmfield(S6,'ssbin'),S6.ssbin==ssbin);
                N(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin) = slength(S7);
                n(looplen-2,looppos,minus2,minus1,plus1,plus2,ssbin,:) = sum(S7.ct,1);
  end,end,end,end,end,end,end,fprintf('\n');
  % SECOND LEVEL OF TABULATION: make list Q of stemstrength series for each looptype+context
  Q=[];qi=0; fs={'looplen','looppos','minus2','minus1','plus1','plus2'}; for fi=1:length(fs),f=fs{fi};Q.(f)=nan(nq,1); end; Q.N=nan(nq,maxssbin); Q.n3=nan(nq,maxssbin,3);
  for looplen=3:maxlooplen, N1=squeeze(N(looplen-2,:,:,:,:,:,:)); n1=squeeze(n(looplen-2,:,:,:,:,:,:,:));
    for looppos=1:looplen, N2=squeeze(N1(looppos,:,:,:,:,:)); n2=squeeze(n1(looppos,:,:,:,:,:,:));
      for minus2=0:4, if minus2==0, N3=squeeze(sum(N2,1)); n3=squeeze(sum(n2,1)); else N3=squeeze(N2(minus2,:,:,:,:)); n3=squeeze(n2(minus2,:,:,:,:,:)); end
        for minus1=0:4, if minus1==0, N4=squeeze(sum(N3,1)); n4=squeeze(sum(n3,1)); else N4=squeeze(N3(minus1,:,:,:)); n4=squeeze(n3(minus1,:,:,:,:)); end
          for plus1=0:4, if plus1==0, N5=squeeze(sum(N4,1)); n5=squeeze(sum(n4,1)); else N5=squeeze(N4(plus1,:,:)); n5=squeeze(n4(plus1,:,:,:)); end
            for plus2=0:4, if plus2==0, N6=squeeze(sum(N5,1)); n6=squeeze(sum(n5,1)); else N6=squeeze(N5(plus2,:)); n6=squeeze(n5(plus2,:,:)); end
              qi=qi+1; for fi=1:length(fs),f=fs{fi};Q.(f)(qi)=eval(f);end; Q.N(qi,:)=N6; Q.n3(qi,:,:)=n6;
  end,end,end,end,end,end, Q=reorder_struct(Q,1:qi); Q.n=sum(Q.n3,3); [Q.rate Q.sd] = ratio_and_sd(Q.n,Q.N); C.Q{ci,1}=Q;
end,toc % 10 min per cohort
% REFORMAT AND SAVE
Z=[]; Z.cohort=rmfield(C,'Q'); Z.stem.ssmin=as_column(ssmin); Z.muttype.name={'C->G';'C->A';'C->T'};
fs={'n3','n','rate','sd'}; Z.loop=rmfield(C.Q{1},fs); for fi=1:length(fs),f=fs{fi};for ci=1:slength(C); Z.loop.(f)(:,:,ci,:)=C.Q{ci}.(f);end,end
save('validation_cohorts_vs_looptypes.tabulated.4.0.mat','Z');

% CURVE FITTING
tmp=load('validation_cohorts_vs_looptypes.tabulated.4.0.mat','Z');X=tmp.Z; % ~5 sec to load
minmut=8; xas=[0:0.5:50 51:1:100 102:2:200]; mhps=[0.30:0.01:0.60];
ssmin = X.stem.ssmin; nss=length(ssmin); x = ssmin; measurepoint=find(x==20);
for ci=1:slength(X.cohort),fprintf('Cohort %d/%d\n',ci,length(X.cohort.name)); npat=X.cohort.npat(ci);
  Q=X.loop; fs={'n3','n','rate','sd'};for fi=1:length(fs),f=fs{fi};Q.(f)=Q.(f)(:,:,ci,:); end;
  W=[]; W.mhp=as_column(mhps); Q.rate=1e6*Q.rate/npat;Q.sd=1e6*Q.sd/npat;
  for wi=1:slength(W),fprintf('%d/%d ',wi,length(W.mhp)); mhp=W.mhp(wi); r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;
    rhp = nan(length(x),length(xas)); for xai=1:length(xas), rhp(:,xai) = r0*2.^(mhp*(x-xas(xai))); end
    for qi=slength(Q):-1:1
      robs = as_column(Q.rate(qi,:)); sdobs = as_column(Q.sd(qi,:)); robs(Q.n(qi,:)<minmut)=nan;
      rinit = robs(1); Q.rinit(qi,1)=rinit; rexp = rinit+rhp; yobs = log10(0.0001+robs); yexp = log10(0.0001+rexp);
      yerr = nansum(bsxfun(@minus,yexp,yobs).^2,1); [err xai] = min(yerr);
      Q.err(qi,1)=err; Q.xa(qi,1)=xas(xai); Q.rexp(qi,:) = rexp(:,xai);
      lastdata = find(~isnan(robs),1,'last'); if isempty(lastdata), lastdata=1; end % cap predictions where data runs out
      Q.rexp(qi,lastdata+1:end)=Q.rexp(qi,lastdata); Q.relrate_exp(qi,:) = Q.rexp(qi,:)/r0; Q.rr = Q.relrate_exp(:,measurepoint);
      Q.rexp(qi,lastdata+1:end)=nan; Q.lastdata(qi,1)=lastdata;
    end,W.err(wi,1)=sum(Q.err); W.Q{wi,1}=Q;
  end,fprintf('\n'); [tmp wi] = min(W.err); X.cohort.mhp(ci,1)=W.mhp(wi); X.cohort.Q{ci,1}=W.Q{wi};
end
save('validation_cohorts_vs_looptypes.tabulated.model.4.0.mat','X');

% Apply model genomewide
tic;load('FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
MODEL=load('validation_cohorts_vs_looptypes.tabulated.model.4.0.mat','X'); MODEL=MODEL.X;
X.site.ss=X.site.nbp+2*X.site.ngc; nssbin=slength(MODEL.stem); X.site.ssbin=min(nssbin,max(1,X.site.ss-3));
X.site.relrate_exp_cohorts = nan(slength(X.site),slength(MODEL.cohort),'single');
for ci=1:slength(MODEL.cohort), fprintf('COHORT %d ',ci); tic
  Q=MODEL.cohort.Q{ci}; Q = reorder_struct(rmfield(Q,{'minus2','plus2'}),Q.minus2==0 & Q.plus2==0);
  rr0=nan(slength(X.site),1,'single'); rr1=rr0;
  for looplen=3:11, i1=find(X.site.looplen==looplen);q1=find(Q.looplen==looplen);
    for looppos=1:looplen, i2=i1(X.site.looppos(i1)==looppos);q2=q1(Q.looppos(q1)==looppos); fprintf('%d/%d ',looppos,looplen);
      for ssbin=1:nssbin, i3=i2(X.site.ssbin(i2)==ssbin);q3=q2(Q.minus1(q2)==0 & Q.plus1(q2)==0); x=ssbin;
        y=q3; rr0(i3)=Q.relrate_exp(y,x);
        for minus1=1:4, i4=i3(X.site.minus1(i3)==minus1);q4=q2(Q.minus1(q2)==minus1);
          for plus1=1:4, i5=i4(X.site.plus1(i4)==plus1);q5=q4(Q.plus1(q4)==plus1);
            y=q5; rr1(i5)=Q.relrate_exp(y,x);
  end,end,end,end,end, fprintf('\n'); toc % ~15 min
  X.site.relrate_exp_cohorts(:,ci) = rr1; idx=find(isnan(rr1)); X.site.relrate_exp_cohorts(idx,ci)=rr0(idx);
end % ~6 min per cohort
save('FINAL_DATASET.with_TpCs.with_validation_models.v4.0.mat','X','-v7.3');
%%%

% how well does model work? (predict rest of patients from just 1 patient)

load('FINAL_DATASET.with_TpCs.with_validation_models.v4.0.mat','X');
pat=X.pat; pat.idx=(1:slength(pat))';
pat = reorder_struct(pat,pat.apobec); % msupe- (<10%) apo+ (>=10%)
pat = sort_struct(pat,'nmut',-1);

X.site.relrate_exp_pat1 = X.site.relrate_exp_cohorts(:,1);  % rate predicted from #1 apobec patient
X.site.relrate_exp_pat2 = X.site.relrate_exp_cohorts(:,2);  % rate predicted from #2 apobec patient
puse=ismember((1:slength(X.pat))',pat.idx(3:end)); muse=puse(X.mut.pat_idx);
X.site.ct_apobec_wgs = histc(X.mut.sidx(muse),1:slength(X.site));  % rate observed at all besides #1+#2

bin=[]; bin.min=[0.2:0.1:2 3:1:10 12 15 20 30 50 70 100]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (X.site.relrate_exp_pat1>=bin.min(i) & X.site.relrate_exp_pat1<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(X.site.ct_apobec_wgs(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(X.site.ct_apobec_wgs); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;
bin1=bin;
bin=[]; bin.min=[0.2:0.1:2 3:1:10 12 15 20 30 50 70 100]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (X.site.relrate_exp_pat2>=bin.min(i) & X.site.relrate_exp_pat2<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(X.site.ct_apobec_wgs(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(X.site.ct_apobec_wgs); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;
bin2=bin;

pr(bin1);pr(bin2);

V={}; V{5}=bin1; V{6}=bin2;
save('validation.tabulated.v4.0.mat','V');
%%%

load('validation.tabulated.v4.0.mat','V');

figure(1)
subplot(1,2,1),cla,hold on
  bin=V{5};
  x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
  scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  xlabel('relative mutation rate (predicted, WGS noncoding, patient #1)','fontsize',20);
  ylabel('relative mutation rate (observed, WGS noncoding, patients #3-178)','fontsize',20);
  line(xlim,xlim);
subplot(1,2,2),cla,hold on
  bin=V{6};
  x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
  scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  xlabel('relative mutation rate (predicted, WGS noncoding, patient #2)','fontsize',20);
  ylabel('relative mutation rate (observed, WGS noncoding, patients #3-178)','fontsize',20);
  line(xlim,xlim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% combine all six validation curves into one file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V={};
tmp=load('validation.tabulated.v1.0.mat','V');V(1:2)=tmp.V(1:2);
tmp=load('validation.tabulated.v2.0.mat','V');V(3)=tmp.V(3);
tmp=load('validation.tabulated.v3.0.mat','V');V(4)=tmp.V(4);
tmp=load('/validation.tabulated.v4.0.mat','V');V(5:6)=tmp.V(5:6);
save('validation.tabulated.1-6.v4.0.mat','V');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIG S10
% six validation curves

% V#                                          __display_order__
% (1) predict WGS coding from WGS noncoding    bottom row left
% (2) predict WXS coding from WGS noncoding    bottom row right
% (3) predict WGS even patients from odd       top row left
% (4) predict WGS even sites from odd          top row right
% (5) predict WGS patients #3-178 from pat #1  middle row left
% (6) predict WGS patients #3-178 from pat #2  middle row right

load('validation.tabulated.1-6.v4.0.mat','V');
figure(1),clf
for i=1:6
  if     i==1, vi=3; xlab='WGS odd patients'; ylab='WGS even patients';
  elseif i==2, vi=4; xlab='WGS odd sites';    ylab='WGS even sites';
  elseif i==3, vi=5; xlab='WGS patient #1';   ylab='WGS patients #3-178';
  elseif i==4, vi=6; xlab='WGS patient #2';   ylab='WGS patients #3-178';
  elseif i==5, vi=1; xlab='WGS noncoding';    ylab='WGS coding';
  elseif i==6, vi=2; xlab='WGS noncoding';    ylab='WXS coding';
  else error('?'); end
  subplot(3,2,i); bin=V{vi};
  rexp=bin.min; robs=bin.relrate; sdobs=bin.relsd;
  x = log2(rexp); y = log2(robs); clr = nansub([0 0 0;1 0 0],1+(x>2|y>2));
  rlo=robs-1.96*sdobs; rhi=robs+1.96*sdobs; idx=find(rlo<0); rlo(idx)=0.1*robs(idx);
  ylo = log2(rlo); yhi = log2(rhi);
  scatter(x,y,50,clr,'filled');ff;set(gca,'fontsize',16);line(xlim,xlim,'color',[0 0 0]);
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  xlabel(xlab,'fontsize',20); ylabel(ylab,'fontsize',20);xlim([-3 8 ]);ylim([-3 8]);
  t=[0.3 1 2 4 10 30 100]; lt=log2(t); set(gca,'xtick',lt,'xticklabel',t,'ytick',lt','yticklabel',t);
  line(xlim,[2 2],'color',[1 0 0]); line([2 2],ylim,'color',[1 0 0]);
  text(2.3,7.2,{'optimal','hairpins'},'color',[1 0 0],'fontsize',16);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TABLES S2-S4

% 63 BASIC LOOP TYPES
load('FINAL.Qs.v3a.1.mat','Qs'); 

Z=[];z=0;
for looplen=3:11, for looppos=1:looplen, Q=Qs{looplen-2}; z=z+1; Z.looplen(z,1)=looplen; Z.looppos(z,1)=looppos;
  for qi=1:slength(Q)-1, Z.(['n' num2str(Q.min(qi))])(z,1) = Q.n(qi,looppos); end
  for qi=1:slength(Q)-1, Z.(['N' num2str(Q.min(qi))])(z,1) = Q.N(qi,looppos); end
  for qi=1:slength(Q)-1, Z.(['rmf' num2str(Q.min(qi))])(z,1) = Q.relrate(qi,looppos); end
end,end
pr(Z)  % --> transfer to Excel sheet

% SEQUENCE CONTEXT

load('FINAL.Qs2.v3a.1.mat','QL3','Q4','Q3','contextL3','context3','context4','Z3','Z4','ZL3','ZL3all');
rc='';rc('ACGT')='TGCA';
for i=1:slength(contextL3),contextL3.name{i}=[rc(contextL3.name{i}(4)) '(' contextL3.name{i}(1:3) ')' contextL3.name{i}(4)];end
for i=1:slength(context3),context3.name{i}=[context3.name{i}(1) '(' context3.name{i}(2:5) ')' rc(context3.name{i}(1))];end
for i=1:slength(context4),context4.name{i}=[rc(context4.name{i}(5)) '(' context4.name{i}(1:4) ')' context4.name{i}(5)];end

% 16 SUBTYPES OF 3/3 LOOPS
Q=QL3;
Z=[];
for z=1:slength(contextL3)
  Z.looptype{z,1}=contextL3.name{z};
  for qi=1:slength(Q), Z.(['n' num2str(Q.min(qi))])(z,1) = Q.n(qi,z); end
  for qi=1:slength(Q), Z.(['N' num2str(Q.min(qi))])(z,1) = Q.N(qi,z); end
  for qi=1:slength(Q), Z.(['rmf' num2str(Q.min(qi))])(z,1) = Q.relrate(qi,z); end
end
pr(Z)

% 64 SUBTYPES OF 3/4 LOOPS
Q=Q3;
Z=[];
for z=1:slength(context3)
  Z.looptype{z,1}=context3.name{z};
  for qi=1:slength(Q), Z.(['n' num2str(Q.min(qi))])(z,1) = Q.n(qi,z); end
  for qi=1:slength(Q), Z.(['N' num2str(Q.min(qi))])(z,1) = Q.N(qi,z); end
  for qi=1:slength(Q), Z.(['rmf' num2str(Q.min(qi))])(z,1) = Q.relrate(qi,z); end
end
Z3=Z;

% 64 SUBTYPES OF 4/4 LOOPS
Q=Q4;
Z=[];
for z=1:slength(context4)
  Z.looptype{z,1}=context4.name{z};
  for qi=1:slength(Q), Z.(['n' num2str(Q.min(qi))])(z,1) = Q.n(qi,z); end
  for qi=1:slength(Q), Z.(['N' num2str(Q.min(qi))])(z,1) = Q.N(qi,z); end
  for qi=1:slength(Q), Z.(['rmf' num2str(Q.min(qi))])(z,1) = Q.relrate(qi,z); end
end
Z4=Z;

Z=concat_structs({Z3,Z4});
pr(Z)




