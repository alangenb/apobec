function R = load_refseq(build)
% R = load_refseq(build)
%
% loads refseq databases from stored files
%
% build can be either a build name, e.g. "hg18", "hg18_v2", "hg19", "mouse" (must not start with "/")
%    OR an absolute path to the Refseq matfile to load (must start with "/")
%
% adjusts all starting coordinates forward by one nucleotide 
% to undo confusing refseq convention of zero-based start, one-based end
%

fprintf('Loading refseq database...\n');

if ~exist('build', 'var') || ~ischar(build)
  fprintf('Assuming hg19.\n');
  build = 'hg19';
 % error('Must provide genome build.  Can be hg17, hg18, hg18_v2, hg19, mm9, OR an absolute path to a Refseq matfile')
end

if strcmpi(build,'hg19gencode')||strcmp(build,'hg19_gencode')
  fprintf('Using hg19 GENCODE.\n');
  build = 'hg19';
  dirname = '/cga/tcga-gsc/home/lawrence/db/hg19/c65e29gencode/';
  matname = [dirname 'R.mat'];
elseif build(1)~='/'
  dirname = ['/xchip/cga/reference/annotation/db/ucsc/' build '/'];
  matname = [dirname 'R.mat'];
else
  if length(build)<5 || ~strcmpi(build(end-3:end),'.mat'), error('should be matfile'); end
  matname = build;
  build = 'hg19';   % process using hg19 codeblock
end

if strcmp(build,'hg19')
  if exist(matname,'file')
  	load(matname,'R')
  else
  	T = read_table([dirname 'refGene.txt'],'%f%s%s%s%f%f%f%f%f%s%s%f%s%s%s%s%s',char(9),0,'whitespace','\b\r');
    [R.id,R.transcript,R.chr,R.strand,R.tx_start,R.tx_end,R.code_start,R.code_end,...
      R.n_exons,R.exon_starts,R.exon_ends,tmp,R.gene,tmp,tmp,R.exon_frames,R.version]  = deal(T.dat{:});
	
	for i=1:slength(R)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
      R.tx_start(i) = R.tx_start(i)+1;
      R.code_start(i) = R.code_start(i)+1;
      R.exon_starts{i} = str2double(split(R.exon_starts{i}(1:end-1),','))+1;
      R.exon_ends{i} = str2double(split(R.exon_ends{i}(1:end-1),','));
      R.exon_frames{i} = str2double(split(R.exon_frames{i}(1:end-1),','));
    end
    fprintf('\n')
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
  	end
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
 	end
    
    save(matname,'R');
  end


elseif strcmp(build,'hg18_v2')
  if exist(matname,'file')
  	load(matname,'R')
  else
  	T = read_table([dirname 'refGene.txt'],'%f%s%s%s%f%f%f%f%f%s%s%f%s%s%s%s%s',char(9),0,'whitespace','\b\r');
    [R.id,R.transcript,R.chr,R.strand,R.tx_start,R.tx_end,R.code_start,R.code_end,...
      R.n_exons,R.exon_starts,R.exon_ends,tmp,R.gene,tmp,tmp,R.exon_frames,R.version]  = deal(T.dat{:});
	
	for i=1:slength(R)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
      R.tx_start(i) = R.tx_start(i)+1;
      R.code_start(i) = R.code_start(i)+1;
      R.exon_starts{i} = str2double(split(R.exon_starts{i}(1:end-1),','))+1;
      R.exon_ends{i} = str2double(split(R.exon_ends{i}(1:end-1),','));
      R.exon_frames{i} = str2double(split(R.exon_frames{i}(1:end-1),','));
    end
    fprintf('\n')
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
  	end
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
 	end
    
    save(matname,'R');
  end
  

elseif strcmp(build,'hg18')

  if exist(matname,'file')
    load(matname,'R');

  else
    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz

    T = read_table([dirname 'refGene.txt'],'%f%s%s%s%f%f%f%f%f%s%s%f%s%s%s%s',...
       char(9),0,'whitespace','\b\r');

    [R.id,R.transcript,R.chr,R.strand,R.tx_start,R.tx_end,R.code_start,R.code_end,...
      R.n_exons,R.exon_starts,R.exon_ends,tmp,R.gene,tmp,tmp,R.exon_frames]  = deal(T.dat{:});

    for i=1:slength(R)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
      R.tx_start(i) = R.tx_start(i)+1;
      R.code_start(i) = R.code_start(i)+1;
      R.exon_starts{i} = str2double(split(R.exon_starts{i}(1:end-1),','))+1;
      R.exon_ends{i} = str2double(split(R.exon_ends{i}(1:end-1),','));
      R.exon_frames{i} = str2double(split(R.exon_frames{i}(1:end-1),','));
    end
    fprintf('\n');

    save(matname,'R');
  end

elseif strcmp(build,'hg17')

  if exist(matname,'file')
    load(matname,'R');

  else

    T = read_table([dirname 'refFlat.txt'],'%s%s%s%s%f%f%f%f%f%s%s',...
      char(9),0,'whitespace','\b\r');

    % wget http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/refFlat.txt.gz

    [R.gene,R.transcript,R.chr,R.strand,R.tx_start,R.tx_end,R.code_start,R.code_end,...
       R.n_exons,R.exon_starts,R.exon_ends]  = deal(T.dat{:});

    for i=1:slength(R)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
      R.tx_start(i) = R.tx_start(i)+1;
      R.code_start(i) = R.code_start(i)+1;
      R.exon_starts{i} = str2double(split(R.exon_starts{i}(1:end-1),','))+1;
      R.exon_ends{i} = str2double(split(R.exon_ends{i}(1:end-1),','));
    end
    fprintf('\n');

    save(matname,'R');
  end
  
elseif strcmp(build,'mm9')
  if exist(matname,'file')
  	load(matname,'R')
  else
  	T = read_table([dirname 'refGene.txt'],'%f%s%s%s%f%f%f%f%f%s%s%f%s%s%s%s%s',char(9),0,'whitespace','\b\r');
    [R.id,R.transcript,R.chr,R.strand,R.tx_start,R.tx_end,R.code_start,R.code_end,...
      R.n_exons,R.exon_starts,R.exon_ends,tmp,R.gene,tmp,tmp,R.exon_frames,R.version]  = deal(T.dat{:});
	
	for i=1:slength(R)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
      R.tx_start(i) = R.tx_start(i)+1;
      R.code_start(i) = R.code_start(i)+1;
      R.exon_starts{i} = str2double(split(R.exon_starts{i}(1:end-1),','))+1;
      R.exon_ends{i} = str2double(split(R.exon_ends{i}(1:end-1),','));
      R.exon_frames{i} = str2double(split(R.exon_frames{i}(1:end-1),','));
    end
    fprintf('\n')
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
  	end
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
 	end
    
    save(matname,'R');
  end
elseif strcmp(build,'canFam2')
  if exist(matname,'file')
  	load(matname,'R')
  else
  	T = read_table([dirname 'refGene.txt'],'%f%s%s%s%f%f%f%f%f%s%s%f%s%s%s%s',char(9),0,'whitespace','\b\r');
    [R.id,R.transcript,R.chr,R.strand,R.tx_start,R.tx_end,R.code_start,R.code_end,...
      R.n_exons,R.exon_starts,R.exon_ends,tmp,R.gene,tmp,tmp,R.exon_frames]  = deal(T.dat{:});
	R.version=zeros(size(R.id));
	for i=1:slength(R)
      if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
      R.tx_start(i) = R.tx_start(i)+1;
      R.code_start(i) = R.code_start(i)+1;
      R.exon_starts{i} = str2double(split(R.exon_starts{i}(1:end-1),','))+1;
      R.exon_ends{i} = str2double(split(R.exon_ends{i}(1:end-1),','));
      R.exon_frames{i} = str2double(split(R.exon_frames{i}(1:end-1),','));
    end
    fprintf('\n')
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
  	end
    
    transcripts_v = cell(slength(R),1);
  	for j = 1:slength(R)
  		transcripts_v{j} = [R.transcript{j} '.' R.version{j}];
 	end
    
    save(matname,'R');
  end  
else
  error('Unknown build %s',build);
end

%% collapse small nucleolar RNA subtypes
R.gene = regexprep(R.gene,'^(SNORD\d*)-\d*$','$1');

%% add txlen, cdlen, protlen
z = zeros(slength(R),1);
R.tx_len = z; R.code_len = z;
for i=1:slength(R)
  for e=1:R.n_exons(i)
    st = R.exon_starts{i}(e);
    en = R.exon_ends{i}(e);
    R.tx_len(i) = R.tx_len(i) + (en-st+1);
    if en<R.code_start(i) || st>R.code_end(i), continue; end
    if st<R.code_start(i), st = R.code_start(i); end
    if en>R.code_end(i), en = R.code_end(i); end
    R.code_len(i) = R.code_len(i) + (en-st+1);
  end
end
R.n_codons = R.code_len / 3;
