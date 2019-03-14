%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FINAL COHORT for revised manuscript
%

%%CHANGE FOLLOWING LINE TO BE FULL PATH WHERE REPO WAS COPIED TO%%
srcpath = '/cga/tcga-gsc/home/alangenb/Lawrence-APOBEC/git/repo/';
addpath(genpath(srcpath));
maxNumCompThreads(8);

pat=makeapn(load_struct('FINAL_DATASET.MIN.v3.2.pat.txt'));
mut=arrayfun(@(x) makeapn(load_struct(['FINAL_DATASET.MIN.v3.2.mut.part' num2str(x) '.txt'])),[1:10]','uni',0);
mut=concat_structs(mut);
ttype=makeapn(load_struct('FINAL_DATASET.MIN.v3.2.ttype.txt'));

X=struct(); X.pat=pat; X.mut=mut; X.ttype=ttype;

save('FINAL_DATASET.MIN.v3.2.mat','X');


%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NMF ON FINAL COHORT
%
% k=8, NMF on rates

% --> (use new file)
%load('FINAL_DATASET.v2.0.mat','X');

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

k=8; randseed=196; X.nmf_rates = perform_nmf(X.nmf_input_rates,k,randseed);

figure(1);clf,display_nmf_legos(X.nmf_rates)  % <------------------- need to manually look at the signatures and figure out which is which



% COULD ADD A STEP HERE:
% --> load the Sanger30 signatures and assign them to our NMFs by cosine similarity-- some should be multiple matches, e.g. APOBEC = Sanger2+Sanger13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANGER signatures 1-30
S = load_struct('/cga/tcga-gsc/home/lawrence/mut/20150107_strand/sanger30.txt');
S = parsein(S,'channel','^(.)\[(.)>(.)\](.)$',{'l','f','t','r'});
bases = {'A','C','G','T'};S.f=listmap(S.f,bases); S.l=listmap(S.l,bases); S.r=listmap(S.r,bases); S.t=listmap(S.t,bases);
idx = find(S.f==4); S.f(idx)=5-S.f(idx); S.t(idx)=5-S.t(idx); tmp=5-S.r(idx); S.r(idx)=5-S.l(idx); S.l(idx)=tmp;
for i=1:slength(S),S.catname{i,1} = [bases{S.f(i)} ' in ' bases{S.l(i)} '_' bases{S.r(i)} ' ->' bases{S.t(i)}]; end
broad = generate_lego_96names(); S.ord = listmap(S.catname,broad); S = sort_struct(S,'ord')
X = {}; for i=1:30, X{i}=str2double(S.(['sig' num2str(i)]));end
%save('/cga/tcga-gsc/home/lawrence/mut/20150107_strand/sanger30.mat','X');
%%%
load('/cga/tcga-gsc/home/lawrence/mut/20150107_strand/sanger30.mat','X');
legos(X,prefix(num2cellstr(1:30),'Signature '));
% saved in /cga/tcga-gsc/home/lawrence/mut/20140218_pancan/Sanger_30_mutational_signatures.pptx
%%%%%%%%%%%%%%




names = {'APOBEC';'UV';'POLE';'MSI';'smoking';'ESO';'aging';'BRCA'};
ord = [7 4 1 3 8 2 5 6];  % <------------------- change this based on what you see above
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

figure(1);clf,display_nmf_legos(X.nmf)

%save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.1.mat','X');
% --> (change name)
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEGO PLOTS

%load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.1.mat','X');
% --> (change name)

figure(1)
ede('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1');
P=[];P.catnames=X.nmf.chan.name;P.add_sanger_plot=true;P.use_lighter_blue=true;
for i=1:slength(X.nmf.factor),disp(X.nmf.factor.name{i});
  clf,lego(X.nmf.chan.nmf(:,i),P);
  set(gcf,'papersize',[3 3],'paperposition',[0.2 0.2 2.6 2.6]);
%  print_to_file(['/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/lego.' X.nmf.factor.name{i} '.pdf']);
% --> (change name)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAP final dataset onto genomic windows

%load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.2.mat','X');
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/hg19_windows.v1.0.mat','W');
X.gene=W.gene; X.win=W.win; clear W
X.mut.zone_idx = get_context(X.mut.chr,X.mut.pos,[srcpath 'ref/zone/']);
X.mut.win_idx = mmw(X.mut,X.win);
puse = (1:slength(X.pat)); muse = find(ismember(X.mut.pat_idx,puse) & X.mut.zone_idx~=3); % nonexon only
X.win.nmut_nonexon_newdata = histc(X.mut.win_idx(muse),1:slength(X.win));
X.win.mutrate_nonexon_newdata = X.win.nmut_nonexon_newdata./(X.win.cov_nonexon*slength(X.pat));
X.win.mutrate_nonexon_newdata(isinf(X.win.mutrate_nonexon_newdata))=nan;  % (smooth behaves differently on nan vs. inf)
%save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET_on_windows.v2.2.mat','X','-v7.3');
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIGURE 1a
% large-scale covariates & mutation rate

%load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET_on_windows.v2.2.mat','X');

chrlen = load_chrlen('hg19'); cen = load_cen('hg19'); window=100000;
chr=17; st=1; en=81e6; ptel=1e5; pcen=1e6; qcen=1e5; qtel=1e6;
Wc = reorder_struct(X.win,X.win.chr==chr);
dat1 = 150+smooth(Wc.rt_extra1,20);
dat2 = 150+max(200,1100-smooth(Wc.expr/500,16));
dat3 = -570+smooth(2.6e8*Wc.mutrate_nonexon_newdata,20);
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
%print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/fig1a.pdf');
%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add more info to mutation object
% --> define patient cohorts

%load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.1.mat','X');
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/gcloud/pats_with_ABdel_classifications.mat','pat');
X.pat = parsein(X.pat,'name','^[^:]+:(.+)$','name0');
X.pat.ABdel = mapacross(X.pat.name0,pat.name,pat.isdel);
X.pat.nmut = histc(X.mut.pat_idx,1:slength(X.pat));
X.pat=order_fields_first(X.pat,{'name','name0','cohort','nmut'});
tic;X.mut = add_llftrr(X.mut,[],[srcpath 'ref/pentamer/']);toc; % ~3min
X.mut = rmfield(X.mut,{'name','num'});
X.pat.frac_apobec = X.pat.nmf_norm(:,1);
X.pat.frac_uv = X.pat.nmf_norm(:,2);
X.pat.frac_pole = X.pat.nmf_norm(:,3);
X.pat.frac_msi = X.pat.nmf_norm(:,4);
X.pat.frac_smoking = X.pat.nmf_norm(:,5);
X.pat.frac_eso = X.pat.nmf_norm(:,6);
X.pat.frac_aging = X.pat.nmf_norm(:,7);
X.pat.frac_brca = X.pat.nmf_norm(:,8);
% Gordenin A3A vs. A3B metric: restricted to Tp(C->K)pA
c2gt = (X.mut.f==2 & X.mut.t~=1);
X.pat.nC = histc(X.mut.pat_idx(c2gt),1:slength(X.pat));
tca2gt = (c2gt & X.mut.l==4 & X.mut.r==1);
X.pat.nRTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==1 | X.mut.ll==3)),1:slength(X.pat));
X.pat.nYTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==2 | X.mut.ll==4)),1:slength(X.pat));
X.pat.RTCA_C = X.pat.nRTCA./X.pat.nC; X.pat.YTCA_C = X.pat.nYTCA./X.pat.nC;
x = X.pat.RTCA_C; y = X.pat.YTCA_C;
X.pat.msupe_neg = (X.pat.frac_msi<0.1&X.pat.frac_smoking<0.1&X.pat.frac_uv<0.1&X.pat.frac_pole<0.1&X.pat.frac_eso<0.1);
X.pat.vanilla = (X.pat.frac_apobec<0.02 & X.pat.msupe_neg);    %  63 APO- MSUPE- patients
X.pat.apobec =  (X.pat.frac_apobec>=0.1 & X.pat.msupe_neg);    % 178 APO+ MSUPE- patients
X.pat.apobec_most_a3a = X.pat.apobec & y>0.281;                %  40 APO+ MSUPE- A3A-most patients
X.pat.apobec_most_a3b = X.pat.apobec & x>0.05 & y<0.116 & (x./(y+0.0315)>=0.655);       %  20 APO+ MSUPE- A3B-most patients
%save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.2.mat','X','-v7.3');
%%%%%%%%%%%%%%%


% output table of patients
%load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.2.mat','X');
pat = X.pat;
pat = sort_struct(pat,{'cohort','ttype','name0'});
pat = keep_fields(pat,{'cohort','name0','ttype','ttype_long','nmut','ABdel','frac_apobec','frac_uv','frac_pole','frac_msi','frac_smoking','frac_eso','frac_aging','frac_brca',...
                    'msupe_neg','vanilla','apobec','RTCA_C','YTCA_C','apobec_most_a3a','apobec_most_a3b'});
pat=rename_field(pat,'name0','name');
%ede('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/supp_data_files');
%save_struct(pat,'/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/supp_data_files/patient_list.txt');
%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A3A vs. A3B PLOTS

%load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.2.mat','X');

% Fig. S4a (Gordenin A3A vs. A3B metric)
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
%print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/fig_S4a.pdf');

% Fig. S4c (lego plots)
%ede('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1');
mut_nonapobec = reorder_struct(X.mut,X.pat.vanilla(X.mut.pat_idx));
mut_apobec = reorder_struct(X.mut,X.pat.apobec(X.mut.pat_idx));
mut_apobec_most_a3b = reorder_struct(X.mut,X.pat.apobec_most_a3b(X.mut.pat_idx));
mut_apobec_most_a3a = reorder_struct(X.mut,X.pat.apobec_most_a3a(X.mut.pat_idx));
figure(1),clf,set(gcf,'papersize',[3 3],'paperposition',[0.2 0.2 2.6 2.6]);
P=[];P.display='rates';P.coverage='genome';P.add_sanger_plot=true;P.use_lighter_blue=true;
clf,legomaf(mut_nonapobec,P);%print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/lego.cohort.non_APOBEC.pdf');
clf,legomaf(mut_apobec,P);%print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/lego.cohort.APOBEC.pdf');
clf,legomaf(mut_apobec_most_a3b,P);%print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/lego.cohort.A3B.pdf');
clf,legomaf(mut_apobec_most_a3a,P);%print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/lego.cohort.A3A.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%
% RUN12 
% --> using final cohort definitions (always exclude other processes from APOBEC cohort)
% --> using process_hairpin_block_v7, which allows specifying gcmin, gcmax as parameters

% compile
mkdir /cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/mcc/phb/v7
cd /cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/mcc/phb/v7
cp /cga/tcga-gsc/home/lawrence/mut/20130429_pancan/mcc/build1/run_matlab.py .
echo "-Xmx16g" > java.opts
reuse .matlab-2016b
mcc -m -C -I /cga/tcga-gsc/home/lawrence/cgal/trunk/matlab/seq -I /cga/tcga-gsc/home/lawrence/cga/trunk/matlab \
    -d /cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/mcc/phb/v7 process_hairpin_block_v7

% run
mutfile = '/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.2.mat';
survdir = '/cga/tcga-gsc/home/lawrence/db/hg19/hairpins/v3a';
gcmin=0.4; gcmax=0.6;
outdir = '/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/wgs/v3/run12'; ede(outdir);
mccdir = '/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/mcc/phb/v7';
mccname = 'process_hairpin_block_v7';
cmds={}; banner='phb';
for blockno=1:32
  cmds{end+1} = ['python run_matlab.py . ' mccname ' ' num2str(blockno) ' ' mutfile ' ' survdir ' ' outdir ' ' num2str(gcmin) ' ' num2str(gcmax)];
end
cmdfile = [outdir '/cmd.txt']; save_lines(cmds,cmdfile); qcmd = ['qsubb ' cmdfile];
qcmd = [qcmd ' --pre=". /broad/software/scripts/useuse;reuse .matlab-2016b;cd ' mccdir '"'];
qcmd = [qcmd ' -cwd -N ' banner ' -j y -l h_rt=2:00:00 -l h_vmem=40g']; qcmd = [qcmd ' -o ' outdir];  % failed with 15g.  jobs sometimes go >1hr
qcmd = [qcmd char(10)]; qcmdfile = [outdir '/qcmd.txt']; save_textfile(qcmd,qcmdfile);

% gather
n=0;N=0;Q=[];nn=[];NN=[]; sss=[0:26]';
indir = '/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/wgs/v3/run12';
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/wgs/v3/run12/stats_all.v1.1.mat','nn','NN','Q','subsets');

prn(ssmin',ssmax',r)

1 0  3     0.886523 0.983310 1.080903 1.005533 0.985504 1.025515 0.946642 0.982103 1.085269 0.892366 0.884584 0.909412
2 4  5     1.113477 1.016690 0.919097 0.994467 1.014496 0.974485 1.053358 1.017897 0.914731 1.107634 1.115416 1.090588
3 6  7     0.863073 0.822852 0.565716 0.822615 0.674710 1.277618 1.147371 0.993938 1.022184 0.870432 0.858924 0.895236
4 8  9     1.096017 0.870437 0.706360 0.886422 0.760412 1.767959 1.139483 1.113654 1.173533 1.082869 1.096589 0.985709
5 10 11    1.063131 0.903889 0.540295 0.902728 0.657036 1.303588 1.130382 1.049367 1.056641 1.053519 1.090669 0.924188
6 12 13    1.493315 0.802368 0.612180 0.802365 0.836542 1.435107 1.112968 1.132892 1.089442 1.466973 1.512316 1.085203
7 14 15    1.946640 0.777216 0.633048 0.775439 0.958659 1.197496 1.010484 1.097268 1.032859 1.808959 2.076498 0.999221
8 16 19    4.044612 0.711367 0.566662 0.770039 0.588241 1.345747 1.173687 0.910245 1.066835 3.859517 4.087215 2.100537
9 20 Inf   7.842426 0.718513 0.516966 0.812600 1.617408 0.988557 1.188874 1.906880 1.237803 7.198451 8.049189 1.605564

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/wgs/v3/run11/stats_all.v1.1.mat','NN','nn','Q','subsets');

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

% FINAL VERSION FOR PAPER (vector-graphics): hairpin-potential decile plots
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/wgs/v3/run11/stats_all.v1.1.mat','NN','nn','Q','subsets');
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
  print_to_file(['/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/hairpins.' subsets{i} '.pdf'],600);
end
pr(subsets,effs,corrs,pvals);

apobec           8.8  0.77  0.015
uv               0.73 -0.90 0.00090
pole             0.48 -0.77 0.014
msi              0.81 -0.83 0.0057
eso              1.6  0.28  0.47
aging            0.96 0.085 0.83
smoking          1.3  0.58  0.099
brca             1.9  0.55  0.12
cohort_nonapobec 1.1  0.50  0.17
cohort_apobec    6.5  0.76  0.017
cohort_A3A       9.0  0.78  0.014
cohort_A3B       1.2  0.54  0.14

% FINAL VERSION FOR PAPER (vector-graphics): replication-timing decile plots
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/wgs/v3/run11/stats_all.v1.1.mat','NN','nn','Q','subsets');
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
  print_to_file(['/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/reptime.' subsets{i} '.pdf'],600);
end
pr(subsets,log2(effs),corrs,pvals);

apobec           -0.079  -0.50 0.14
uv               0.64    0.99  1.81e-08
pole             1.1     0.99  1.75e-08
msi              0.69    0.99  5.33e-09
eso              4.4     0.95  2.89e-05
aging            1.1     0.97  2.89e-06
smoking          2.7     0.96  1.23e-05
brca             0.74    0.99  3.74e-08
cohort_nonapobec 0.64    0.94  4.70e-05
cohort_apobec    0.33    0.93  0.00011
cohort_A3A       -0.16   -0.85 0.0019
cohort_A3B       -0.0043 0.42  0.23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MAP MUTATIONS TO LIST OF TpCs

% hairpins list v3b (looppos = 1 - looplen)
tpcfile = '/cga/tcga-gsc/home/lawrence/db/hg19/hairpins/v3a/all.TpCs_only.mat';
mutfile = '/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.2.mat';
tic;load(tpcfile,'X');toc % 40 sec
tic;tmp=load(mutfile);toc % 10 sec
fs={'pat','ttype','nmf','mut'};for fi=1:length(fs),f=fs{fi};X.(f)=tmp.X.(f);end
tic;X.mut.sidx = multimap(X.mut,X.site,{'chr','pos'});toc % 12 min
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X','-v7.3');

% hairpins list v3b (looppos = 0 - looplen+1)
tpcfile = '/cga/tcga-gsc/home/lawrence/db/hg19/hairpins/v3b/all.TpCs_only.mat';
mutfile = '/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.v2.2.mat';
tic;load(tpcfile,'X');toc % 40 sec
tic;tmp=load(mutfile);toc % 10 sec
fs={'pat','ttype','nmf','mut'};for fi=1:length(fs),f=fs{fi};X.(f)=tmp.X.(f);end
tic;X.mut.sidx = multimap(X.mut,X.site,{'chr','pos'});toc % 12 min
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3b.v2.2.mat','X','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figs.3b,S6

% TABULATE 

tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec to load --> 9.3Gb in mem
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.Qs.v3a.1.mat','Qs');

% PLOT for Fig. 3b
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.Qs.v3a.1.mat','Qs');
figure(4),clf,looplen=4;Q=Qs{looplen-2};Q=reorder_struct(Q,1:slength(Q)-1);rate=Q.relrate;sd=Q.relsd;ratehi=rate+1.96*sd;ratelo=rate-1.96*sd;
b=barweb(rate,ratelo,ratehi);ff;colormap(parula);set(gca,'xtick',1:slength(Q),'xticklabel',Q.label,'fontsize',14,'linewidth',1.5);
for bi=1:length(b.bars),set(b.bars(bi),'linewidth',1);end;xlim([0.3 slength(Q)+0.7]);line(xlim,[1 1],'linestyle',':','color',[0 0 0]);
set(gcf,'papersize',[10 7.5],'paperposition',[0.2 0.2 10-0.4 7.5-0.4]);
print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/fig.3b.pdf');

% PLOT for Fig. S6
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.Qs.v3a.1.mat','Qs');
figure(5),clf,looplen=3;Q=Qs{looplen-2};Q=reorder_struct(Q,1:slength(Q)-1);rate=Q.relrate;sd=Q.relsd;ratehi=rate+1.96*sd;ratelo=rate-1.96*sd;
b=barweb(rate,ratelo,ratehi);ff;colormap(parula);set(gca,'xtick',1:slength(Q),'xticklabel',Q.label,'fontsize',14,'linewidth',1.5);
for bi=1:length(b.bars),set(b.bars(bi),'linewidth',1);end;xlim([0.3 slength(Q)+0.7]);line(xlim,[1 1],'linestyle',':','color',[0 0 0]);
set(gcf,'papersize',[10 7.5],'paperposition',[0.2 0.2 10-0.4 7.5-0.4]);
print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/fig.S6.pdf');

% PLOT for Fig. 3d
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.Qs.v3a.1.mat','Qs');
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
print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/fig.3d.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figs.3f,S8

% TABULATE 

tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec to load --> 9.3Gb in mem
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.Qs2.v3a.1.mat','QL3','Q4','Q3','contextL3','context3','context4','Z3','Z4','ZL3','ZL3all');

%  check numbers
pr(Z4);pr(Z3);pr(ZL3);pr(ZL3all)
%%%

% PLOT  OLD Fig. 3f, S8

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.Qs2.v3a.1.mat','Z3','Z4','ZL3','ZL3all');

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
% --> LOOKS PERFECT!


% PLOT  NEW Fig. 3f

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.Qs2.v3a.1.mat','QL3','contextL3');
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
print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/fig.3f.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SUBMISSION due Weds 2019-03-06

% 1. Assemble new supp info
% 2. Revise text
% 3. Revise rebuttal
% 4. Assemble sourcecode + files to share
% 5. Assemble Alexandrov/Nik-Zainal-only dataset to distribute with code.
%    --> See how well it replicates the various figures.


% FIGURE 1
%   Fig. 1a: update red curve
%            update bottom plots
%            remove A-F letters

% SUPP INFO
% --> make sure all of Remi's new data gets into the supp info (add them into the powerpoint)

% FURTHER ANALYSES REQUESTED BY REVIEWERS




% QUANTITATIVE MODEL
% --> simplest path to a new Fig.4 is as follows:
%   % 3. Define cohorts in advance, once, for all analyses:
%    --> one cohort per NMF factor, containing all patients that are >=50% that factor.  for non-APOBEC cohorts, require frac_APOBEC<
%        require patients to be >50% APOBEC and <10% for each of POLE/MSI/UV/SMOKING/ESO (AGING/BRCA are ok)
%    --> APOBEC-  means APOBEC<0.05
%        APOBEC+  means APOBEC>=0.10
%        APOBEC++ means APOBEC>=0.50
%        "APOBEC" cohort means "APOBEC+ POLE- MSI- UV- SMOKING- ESO-" for the purposes of 


% MAJOR TO-DO
% 1. Replace all plots in Figs.1,S1,S4,S6,S8
%    --> possibly combine S6 and S8 into a single figure.  This frees up S8 for other stuff Remi needs to add for the rebuttal


% 4. revise the A3A vs. A3B analysis: require patients to be >50% APOBEC and <10% for each of POLE/MSI/UV/SMOKING/ESO (i.e. at least  AGING+BRCA 
% 4. Quantitative model: for now, commit to using just (looplen,looppos,minus1,plus1,stemstrength), i.e. ignore minus2,plus2
%    -->

% 4. Show that non-APOBEC (green dots) hairpin flatness remains when analyzing just the Tp(C->G)


% --> add a generalized 3/3 hairpin to the whitespace in fig. 3f

% CODE AVAILABILITY
% Laura:
%   Code must be provided and deposited in a public repository that allows for freezing the code (IE Zenodo, not Github) and declared in the acknowledgments.
% Checklist:
%   Computer code stored in GitHub should be permanently archived in a repository such as Zenodo or Code Ocean
%   and cited in the reference list with an associated DOI. 
%   For guidance, see https://help.github.com/articles/referencing-and-citing-content/.


% MINOR TO-DO

% 1. Slightly increase size of bars in Fig. 3d (reduce x-axis gaps)
% 2. Replace "smoking/aging" with "SMOKING/AGING" everywhere in the text
% 3. Increase font size, decrease stroke size in all hairpin diagrams -- to make more readable

% FUN MINOR TO-DO (save these for after Weds.)

% 1. Compare stem strength metrics: which best correlates to the in vitro data?  (project for new person in lab)
% 2. Gordenin's A3A vs. A3B metric: is currently restricted to Tp(C->K)pA.  Does it improve with inclusion of Tp(C->K)pT? Tp(C->G)pG?
% 3. Kataegis: focus on the top handful of most-kataegic samples.  
%         --> do the APOBEC clusters include many non-TpC sites?  
%         --> if so, this directly proves that APOBEC can stimulate non-TpC mutations!
% 4. Iteratively test each nested model: e.g. does minus2 ever matter by itself? conditionally dependent on minus1?
% 5. Estimate a per-sample parameter, fmutA3A = what fraction of the total mutational burden is due to A3A?
%         --> probably can do this from the slope of 3/3.  (test whether including larger loops ever helps)
% 6. Look at nonmethylated CpGs (in our favor: APOBEC seems to be uniquely uniform across large-scale genomic covariates)
% 7. Currently we consider each (looplen,looppos) separately. Look at assigning each a scalar multiplier that can be calculated analytically from (looplen,looppos)
%         --> first, remap to (looplen,tpos),  then estimate a global tpos_optimal
%         --> imagine the response curve as a gaussian over the tpos0,1 interval
%         --> plot: one curve per looplen (use colors), with  y=rrmax (ss>=19 in 3/3) vs. x=(looppos-1)/looplen;
% 8. Does it matter whether we calculate frac_apobec using NMF on (rates or counts) ?
% 9. Expand to GC=30%-70%, GC=20%-80%, GC=10%-90%; make an animation of these expanding subsets.  Also look at each decile individually.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WXS hotspots
% from MC3
%
% Need to update the list of driver genes... download the latest Sanger Cosmic CGC list
%   from https://cancer.sanger.ac.uk/census

% this was the link we downloaded from:
https://cancer.sanger.ac.uk/cosmic/census/all?home=y&name=all&tier=&sEcho=1&iColumns=20&sColumns=&iDisplayStart=0&iDisplayLength=25&mDataProp_0=0&sSearch_0=&bRegex_0=false&bSearchable_0=true&bSortable_0=true&mDataProp_1=1&sSearch_1=&bRegex_1=false&bSearchable_1=true&bSortable_1=true&mDataProp_2=2&sSearch_2=&bRegex_2=false&bSearchable_2=true&bSortable_2=true&mDataProp_3=3&sSearch_3=&bRegex_3=false&bSearchable_3=true&bSortable_3=true&mDataProp_4=4&sSearch_4=&bRegex_4=false&bSearchable_4=true&bSortable_4=true&mDataProp_5=5&sSearch_5=&bRegex_5=false&bSearchable_5=true&bSortable_5=true&mDataProp_6=6&sSearch_6=&bRegex_6=false&bSearchable_6=true&bSortable_6=true&mDataProp_7=7&sSearch_7=&bRegex_7=false&bSearchable_7=true&bSortable_7=true&mDataProp_8=8&sSearch_8=&bRegex_8=false&bSearchable_8=true&bSortable_8=true&mDataProp_9=9&sSearch_9=&bRegex_9=false&bSearchable_9=true&bSortable_9=true&mDataProp_10=10&sSearch_10=&bRegex_10=false&bSearchable_10=true&bSortable_10=true&mDataProp_11=11&sSearch_11=&bRegex_11=false&bSearchable_11=true&bSortable_11=true&mDataProp_12=12&sSearch_12=&bRegex_12=false&bSearchable_12=true&bSortable_12=true&mDataProp_13=13&sSearch_13=&bRegex_13=false&bSearchable_13=true&bSortable_13=true&mDataProp_14=14&sSearch_14=&bRegex_14=false&bSearchable_14=true&bSortable_14=true&mDataProp_15=15&sSearch_15=&bRegex_15=false&bSearchable_15=true&bSortable_15=true&mDataProp_16=16&sSearch_16=&bRegex_16=false&bSearchable_16=true&bSortable_16=true&mDataProp_17=17&sSearch_17=&bRegex_17=false&bSearchable_17=true&bSortable_17=true&mDataProp_18=18&sSearch_18=&bRegex_18=false&bSearchable_18=true&bSortable_18=true&mDataProp_19=19&sSearch_19=&bRegex_19=false&bSearchable_19=true&bSortable_19=true&sSearch=&bRegex=false&iSortCol_0=0&sSortDir_0=asc&iSortingCols=1&export=tsv

% this link also works:
https://cancer.sanger.ac.uk/cosmic/census/all?home=y&name=all&export=tsv

%   original citation was: https://www.nature.com/articles/nrc1299
% --> downloaded Census_allTue Mar  5 12_59_32 2019.tsv
% to: /cga/tcga-gsc/home/lawrence/db/cgc/Census_allTue_Mar_5_12_59_32_2019.tsv
%
% About the CGC list, Laura Zahn says:
% Please make this a numbered reference that gives where this information can be obtained and use the number here.
%

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/all_12086_WXS_hotspots.2.3.mat','S');
f=grep('relrate|qidx|column',fieldnames(S));S=rename_fields(S,f,prefix(f,'old_'));

% annotate with status in latest CGC driver list
C=mals('/cga/tcga-gsc/home/lawrence/db/cgc/Census_allTue_Mar_5_12_59_32_2019.tsv');
C=reorder_struct(C,grepm('Mis|N|S',C.MutationTypes)); % --> 348 mutation-driven genes
S=rmfield(S,'driver'); S.cgc = ismember(S.gene,C.GeneSymbol);

save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/all_12086_WXS_hotspots.2.4.mat','S');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FINAL QUANTITATIVE MODELING

tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
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
C.name{end+1,1}='APOBEC>=90%'; C.puse(end+1,:)=ismember((1:slength(X.pat)),pat.idx(pat.frac_apobec>=0.9));
C.name{end+1,1}='APOBEC>=50%'; C.puse(end+1,:)=ismember((1:slength(X.pat)),pat.idx(pat.frac_apobec>=0.5));
C.name{end+1,1}='APOBEC>=10%'; C.puse(end+1,:)=ismember((1:slength(X.pat)),pat.idx(pat.frac_apobec>=0.1));
C.name{end+1,1}='APOBEC>=3%'; C.puse(end+1,:)=ismember((1:slength(X.pat)),pat.idx(pat.frac_apobec>=0.03));
C.npat = sum(C.puse,2);
for ci=1:length(C.name),fprintf('%d/%d ',ci,length(C.name));
  C.muse{ci,1} = C.puse(ci,X.mut.pat_idx)==1; C.nmut(ci,1) = sum(C.muse{ci});
end,fprintf('\n');pr(C)

name_______ puse  npat muse  nmut
APOBEC>=90% {...} 5    {...} 184310
APOBEC>=50% {...} 76   {...} 1110525
APOBEC>=10% {...} 178  {...} 1399977
APOBEC>=3%  {...} 244  {...} 1478964

% TABULATE 
tic
maxlooplen=11; maxssbin=20; nq = round(1.1*sum(3:maxlooplen)*(4.^4));
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_cohorts_vs_looptypes.tabulated.1.0.mat','Z');
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CURVE FITTING FOR EACH SERIES

tmp=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_cohorts_vs_looptypes.tabulated.1.0.mat','Z');X=tmp.Z; % ~5 sec to load
minmut=8; xas=[0:0.5:50 51:1:100 102:2:200]; mhps=0.45;
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SHOW FITS
% on triangular arrangement of subplots

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X');
ci=4; Q=X.cohort.Q{ci}; npat=X.cohort.npat(ci);r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% APPLY MODEL TO WXS HOTSPOTS TO DISTINGUISH
% PASSENGERS VS. DRIVERS
%
% --> annotate the list of MC3 hotspots with relrate_exp from the model.
% --> If there's enough data, use plus1/minus1 to refine relrate_exp.  Otherwise just use looplen/looppos.

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/all_12086_WXS_hotspots.2.4.mat','S');

MODEL=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X');MODEL=MODEL.X;
ci=3; Q=MODEL.cohort.Q{ci}; Q = reorder_struct(rmfield(Q,{'minus2','plus2'}),Q.minus2==0 & Q.plus2==0);

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

save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.all_12086_WXS_hotspots.3.0.mat','S');
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% update Table S1
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.all_12086_WXS_hotspots.3.0.mat','S');
save_struct(keep_fields(S,{'gene','cgc','relrate_exp'}),'/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.all_12086_WXS_hotspots.3.0.for_TableS1.txt');
% --> integrated CGC membership and relrate_exp into existing table
%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIGURE 4

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL.all_12086_WXS_hotspots.3.0.mat','S');
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
S.showname(grepm('HIST',S.gene))=0;P=[];P.xpad=1; P.ypad=1;
idx=find(S.showname); textfit(x(idx)+(0.010+0.004*tested(idx))*diff(xlim),y(idx),S.gene(idx),'fontsize',fsz(idx),'fontcolor',clr(idx,:),P);
% --> NOTE: need to manually enlarge figure window and then repeat rendering, to get textit spacing right, *then* write to PDF.
set(gcf,'papersize',[12 10],'paperposition',[0.2 0.2 12-0.4 10-0.4]);
print_to_file('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/finalfigs/sampledata/v1/fig.4.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% APPLY MODEL GENOMEWIDE
%
% Annotate each TpC site withrelrate_exp
% --> in cases where there is enough data to use minus1,plus1, use those to refine rate.
% --> otherwise, just use looplen/looppos

tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec

MODEL=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_cohorts_vs_looptypes.tabulated.model.1.2.mat','X'); MODEL=MODEL.X;
ci=3; Q=MODEL.cohort.Q{ci}; Q = reorder_struct(rmfield(Q,{'minus2','plus2'}),Q.minus2==0 & Q.plus2==0);

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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs.with_model.v1.0.mat','X','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% List of top A3A hairpin sites for Adam's chrfig plots
%
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs.with_model.v1.0.mat','X');
S=X.site;
puse=(X.pat.apobec); muse=puse(X.mut.pat_idx); % msupe- apo+
S.ct_apobec_wgs = histc(X.mut.sidx(muse),1:slength(S));
S.idx_orig = (1:slength(S))'; S=off(S,'idx_orig');
S = reorder_struct(S,S.relrate_exp>=10);
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.top_hairpins.v1.0.mat','S');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%% CURVE FITTING CODE %%%%%%%%%%%%%%
tmp=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_cohorts_vs_looptypes.tabulated.1.0.mat','Z');X=tmp.Z; % ~5 sec to load
minmut=8; xas=[0:0.5:50 51:1:100 102:2:200]; mhps=0.45;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% VALIDATION OF MODELING

% (1) Predict WGS coding rates from WGS noncoding model
% (2) Predict WXS coding rates from WGS noncoding model

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs.with_model.v1.0.mat','X');
puse=(X.pat.apobec); muse=puse(X.mut.pat_idx); % msupe- apo+
Sx = X.site; Sx.ct_apobec_wgs = histc(X.mut.sidx(muse),1:slength(Sx));
Sx.idx_orig = (1:slength(Sx))'; Sx=off(Sx,'idx_orig'); Sx = reorder_struct(rmfield(Sx,'zone'),Sx.zone==3);
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.exonic_TpCs.with_model.v1.0.mat','Sx');
%%%
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.exonic_TpCs.with_model.v1.0.mat','Sx');
load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/MC3.exome_TpC_sites.spliceflank10.1.4.Gflipped.compact.1.0.mat','M');
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

save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.v1.0.mat','Sx');
%%%

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.v1.0.mat','Sx');

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

save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v1.0.mat','V');

% more conservative version: go out to just 50
V={};

bin=[]; bin.min=[0.2:0.1:2 3:1:10 12 15 20 30 50]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (Sx.relrate_exp>=bin.min(i) & Sx.relrate_exp<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(Sx.ct_apobec_wgs(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(Sx.ct_apobec_wgs); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;

V{1}=bin;

bin=[]; bin.min=[0.2:0.1:3 3.5:0.5:10 12:1:20 25 30 50]'; bin.max=[bin.min(2:end);inf];
for i=1:slength(bin)
  ii = (Sx.relrate_exp>=bin.min(i) & Sx.relrate_exp<bin.max(i));
  bin.N(i,1)=sum(ii); bin.n(i,1)=sum(Sx.ct_apobec_wxs(ii));
end, [bin.rate bin.sd]  = ratio_and_sd(bin.n,bin.N);
r0 = mean(Sx.ct_apobec_wxs); bin.relrate = bin.rate/r0; bin.relsd = bin.sd/r0;

V{2}=bin;

save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated_to_50x.v1.0.mat','V');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next:
% (3) split patients in half, predict one half from the other

% REPEAT TABULATION AND CURVE FITTING FOR THESE

tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
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
X.mut.odd = mod(X.mut.pos,2)==1;
X.mut.even = mod(X.mut.pos,2)==0;

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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.2.0.mat','Z');

% CURVE FITTING
tmp=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.2.0.mat','Z');X=tmp.Z; % ~5 sec to load
minmut=8; xas=[0:0.5:50 51:1:100 102:2:200]; mhps=[0.30:0.01:0.60]; %mhps=0.45;
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.model.2.0.mat','X');

pr(X.cohort)
%%%

% Apply models genomewide
tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
MODEL=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.model.2.0.mat','X'); MODEL=MODEL.X;
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs.with_validation_models.v2.0.mat','X','-v7.3');
%%%

% how well does model work? (predict even patients from odd patients)

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs.with_validation_models.v2.0.mat','X');
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v2.0.mat','V');
%%%

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v2.0.mat','V');

figure(1)
bin=V{3};
x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff
ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
xlabel('relative mutation rate (predicted, WGS noncoding, odd patients)','fontsize',20);
ylabel('relative mutation rate (observed, WGS noncoding, even patients)','fontsize',20);
line(xlim,xlim);
% --> looks perfect!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NEXT: do the even/odd sites test properly 
% (did even/odd mutations last time, wasn't thinking clearly)

% (4) predict even sites from odd sites

% REPEAT TABULATION AND CURVE FITTING

tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.3.0.mat','Z');

% CURVE FITTING
tmp=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.3.0.mat','Z');X=tmp.Z; % ~5 sec to load
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.model.3.0.mat','X');

pr(X.cohort)

% Apply model genomewide
tic;load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs_v3a.v2.2.mat','X');toc   % 50 sec
MODEL=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation_cohorts_vs_looptypes.tabulated.model.3.0.mat','X'); MODEL=MODEL.X;
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs.with_validation_models.v3.0.mat','X','-v7.3');
%%%

% how well does model work? (predict even sites from odd sites)

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/FINAL_DATASET.with_TpCs.with_validation_models.v3.0.mat','X');
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
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v3.0.mat','V');
%%%

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v3.0.mat','V');

figure(1)
bin=V{4};
x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff
ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
xlabel('relative mutation rate (predicted, WGS noncoding, odd sites)','fontsize',20);
ylabel('relative mutation rate (observed, WGS noncoding, even sites)','fontsize',20);
line(xlim,xlim);
% --> looks great!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combine all four validation curves into one file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V={};
tmp=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v1.0.mat','V');V(1:2)=tmp.V(1:2);
tmp=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v2.0.mat','V');V(3)=tmp.V(3);
tmp=load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.v3.0.mat','V');V(4)=tmp.V(4);
save('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.1-4.v4.0.mat','V');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make combined version for new Fig.S8 that has all four validation curves

load('/cga/tcga-gsc/home/lawrence/apobec/20170117_rnaed/validation.tabulated.1-4.v4.0.mat','V');

figure(1),clf
subplot(2,2,1),bin=V{1};
  x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
  scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff;set(gca,'fontsize',16);
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  xlabel('predicted rate, WGS noncoding','fontsize',20);
  ylabel('observed rate, WGS coding','fontsize',20);
  line(xlim,xlim,'color',[0 0 0]);ylim([0 100]);
subplot(2,2,2),bin = V{2};
  x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
  scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff;set(gca,'fontsize',16);
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  xlabel('predicted rate, WGS noncoding','fontsize',20);
  ylabel('observed rate, WXS coding','fontsize',20);
  line(xlim,xlim,'color',[0 0 0]);ylim([0 180]);
subplot(2,2,3),bin=V{3};
  x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
  scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff;set(gca,'fontsize',16);
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  xlabel('predicted rate, WGS odd patients','fontsize',20);
  ylabel('observed rate, WGS even patients','fontsize',20);
  line(xlim,xlim,'color',[0 0 0]);
subplot(2,2,4),bin=V{4};
  x = bin.min; y = bin.relrate; ylo = bin.relrate-1.96*bin.relsd; yhi = bin.relrate+1.96*bin.relsd;
  scatter(bin.min,bin.relrate,50,[0 0 0],'filled');ff;set(gca,'fontsize',16);
  ebw=0.2;for j=1:length(x),line([1 1]*x(j),[ylo(j) yhi(j)],'color',[0 0 0],'linewidth',1);line(x(j)+[-1 1]*ebw,[1;1]*[yhi(j) ylo(j)],'color',[0 0 0],'linewidth',1);end
  xlabel('predicted rate, WGS odd sites','fontsize',20);
  ylabel('observed rate, WGS even sites','fontsize',20);
  line(xlim,xlim,'color',[0 0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


