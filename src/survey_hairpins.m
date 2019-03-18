function survey_hairpins(outdir,blockno)

if nargin~=2, error('should take two input arguments: outdir, blockno'); end

if ~ischar(outdir), error('outdir should be a string'); end
if ischar(blockno), blockno=str2double(blockno); end

ede(outdir);

%spliceflank=0;   % ORIGINALLY WAS THIS
spliceflank=20;   % SHOULD BE THIS (tweaked this manually afterwards)
X=[]; tmp=load('/cga/tcga-gsc/home/lawrence/db/hg19/hg19_genome_blocks.v1.0.mat','B'); X.block=tmp.B; X.gene = load_genes('hg19gencode',spliceflank);
if blockno<1 || blockno>slength(X.block) || blockno~=round(blockno), error('invalid blockno'); end

X.loop=[];
X.loop.len = uint8([repmat(3,1,3) repmat(4,1,4) repmat(5,1,5) repmat(6,1,6) repmat(7,1,7) repmat(8,1,8) repmat(9,1,9) repmat(10,1,10) repmat(11,1,11)]');
X.loop.pos = uint8([1:3 1:4 1:5 1:6 1:7 1:8 1:9 1:10 1:11]');
X.loop.pos_flip = X.loop.len+1-X.loop.pos; X.loop.flip_idx = multimap(X.loop,X.loop,{'len','pos'},{'len','pos_flip'});
X.loop.tpos = (double(X.loop.pos)-1)-(double((X.loop.len+1))/2); X.loop.good = (abs(X.loop.tpos)<=0.5);
X.loop.pos_double = double(X.loop.pos); X.loop.len_double = double(X.loop.len); maxstem=250; window=maxstem+max(X.loop.len_double);

outfile = [outdir '/block' num2str(blockno) '.mat'];
write_textfile('IN_PROGRESS',outfile);
fprintf('BLOCK %d\n',blockno);

X.block.this = (1:slength(X.block))'==blockno; chr=X.block.chr(blockno); st=X.block.st(blockno); en=X.block.en(blockno); len=en-st+1;
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
X.site.nbp = zeros(len,slength(X.loop),'uint8'); X.site.ngc = X.site.nbp; fprintf('Mb:');
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




