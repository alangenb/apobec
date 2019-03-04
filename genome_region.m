function r = genome_region(c, s, e, P)
% genome_region(chromosome, start, end, build)
%
% P.build can be hg17, hg18, hg19, mm9, OR absolute path to directory of chr*.txt files
% P.refir is absolute path to directory of <build>_info.txt and chr*.txt files

if nargin==2
  e=s;
end

if ischar(s), s=str2double(s); end
if ischar(e), e=str2double(e); end

if ~exist('P','var'), P=[]; end

if ischar(P)
    build=P;
    P=[];
    P.build=build;
end

if ~exist('build','var')
    build = [];
end

if strcmp(build,'hg19gencode'), build='hg19'; end


P = impose_default_value(P,'refdir',[]);



if ~exist('build', 'var') || ~ischar(build)
  fprintf('Build not provided... Using hg19 by default\n');
  build = 'hg19';
%    error('Must provide genome build.  Can be hg17, hg18, hg19, mm9, OR absolute path to directory of chr*.txt files')
end

% METHODS AVAILABLE:
% if ReferenceInfoObj has been initialized, use it.
% if ReferenceInfoObj has not been initialized, use legacy methods.

try

  maxnum = ReferenceInfoObj.getMaxNum(build);

  % OK, it worked: we can use ReferenceInfoObj

  if ((~isempty(P.build)) &  (~isempty(P.refdir) ) )
    
    % CS 7/2012
    ReferenceInfoObj.init(P.refdir)
    c=char(ReferenceInfoObj.num2chrom(c,P.build));
    build=P.build;
    dirname=P.refdir;
    
  elseif (~isempty(P.build) )
    
    % CS 7/2012
    r1=char(ReferenceInfoObj.getReferenceFile(c));
    build=P.build;
    [p1, f1, x1] = fileparts(r1); 
    dirname=[p1 '/'];
    
  end

catch me

  % ReferenceInfoObj not available: use legacy methods
    
  if build(1) ~= '/'
    if strcmp(build,'hg17')
      dirname = ['/xchip/tcga/gbm/analysis/lawrence/genome/' build '/'];
    else
      dirname = ['/xchip/cga/reference/annotation/db/ucsc/' build '/'];
    end
  else
    dirname = build;
    if dirname(end)~='/', dirname = [dirname '/']; end
    %  clear build;
  end
    
end



if isnumeric(c)
    if length(c)>1, error('multiple chromosomes not supported'); end
    c = chrlabel(c,build);

%    if strcmp(build,'mm9')||strcmp(build(end-3:end),'mm9/')||strcmp(build(end-2:end),'mm9')
%        if c==20, c = 'X';
%        elseif c==21, c = 'Y';
%        elseif c==22, c = 'M'; end
%    end
%    if strcmp(build,'canFam2')||strcmp(build(end-7:end),'canFam2/')||strcmp(build(end-6:end),'canFam2')
%        if c==20, c = 'X';
%        elseif c==21, c = 'Y';
%        elseif c==22, c = 'M'; end
%    end
%    if c==23, c = 'X';
%    elseif c==24, c = 'Y';
%    elseif c==25, c = 'M';
%    else c = num2str(c); end
end

if iscell(c) && length(c)==1, c=c{1}; end
if ~ischar(c), error('chromsome should be numeric or a string'); end
if ~strncmp(c,'chr',3), c = ['chr' c]; end

%if strcmp(build,'mm9')||strcmp(build(end-3:end),'mm9/')||strcmp(build(end-2:end),'mm9')
%    if strcmp(c,'chr20'), c = 'chrX'; end
%    if strcmp(c,'chr21'), c = 'chrY'; end
%end
%if strcmp(c,'chr23'), c = 'chrX'; end
%if strcmp(c,'chr24'), c = 'chrY'; end
%ok = {'chr0';'chr10_random';'chr10';'chr11_random';'chr11';'chr12';'chr13_random'; ...
%    'chr13';'chr14';'chr15_random';'chr15';'chr16_random';'chr16';'chr17_random'; ...
%    'chr17';'chr18_random';'chr18';'chr19_random';'chr19';'chr1_random';'chr1'; ...
%    'chr20';'chr21_random';'chr21';'chr22_h2_hap1';'chr22_random';'chr22';'chr23'; ...
%    'chr24';'chr2_random';'chr2';'chr3_random';'chr3';'chr4_random';'chr4'; ...
%    'chr5_h2_hap1';'chr5_random';'chr5';'chr6_cox_hap1';'chr6_qbl_hap2';'chr6_random'; ...
%    'chr6';'chr7_random';'chr7';'chr8_random';'chr8';'chr9_random';'chr9';'chrM'; ...
%    'chrX_random';'chrX';'chrY'};
%
%if ~ismember(c,ok), error('%s.txt is not a valid genome file',c); end

fname = [dirname c '.txt'];

if ~exist('e','var'), e=s; end
if any(size(s)>1) | any(size(e)>1)   % vector mode
    if length(s)~=length(s(:)) | length(e)~=length(e(:)), error('"s" and "e" can be vectors but not matrices'); end
    if length(s)~=length(e), error('"s" and "e" must be same length'); end
end
sv = s;
ev = e;
r = cell(length(sv),1);

% find file
d = dir(fname);
if length(d)~=1, error('Couldn''t find genome file %s',fname); end
filesize = d.bytes;

% open file
for attempt=1:10
    try
        f = fopen(fname, 'rt');
        if (f==-1), error(['Could not open ' fname]); end
        break
    catch me
        fprintf('Failed to open %s\n',fname);
        disp(me); disp(me.message);
        fprintf('Waiting ten seconds and trying again...\n');
        pause(10);   % wait ten seconds and try again
    end
end   % next try

% look up sequence(s)
for i=1:length(sv)
    s = sv(i);
    e = ev(i);

    if e<s, error('end<start'); end
    if s<1, s = 1; fprintf('WARNING: trimming <start> to 1\n'); end
    if e<1, e = 1; fprintf('WARNING: trimming <end> to 1\n'); end
    if s>filesize, s = filesize; fprintf('WARNING: trimming <start> to <filesize>\n'); end
    if e>filesize, e = filesize; fprintf('WARNING: trimming <end> to <filesize>\n'); end
    
    status = fseek(f, s-1, 'bof');
    if (status ~= 0), error(['Could not seek to ' num2str(s) '-1 in ' fname]); end
    sz = e-s+1;
    r{i} = fgets(f,sz);
end % next i

fclose(f);

if length(sv)==1, r = r{1}; end   % cell->char for single queries
    
    
    
