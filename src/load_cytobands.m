function C = load_cytobands(build)

if ~exist('build','var')
  fprintf('Assuming hg18');
  build = 'hg18';
end

if strcmpi(build,'hg17')
  dr = '/xchip/cga/reference/annotation/db/ucsc/hg17';
elseif strcmpi(build,'hg18')
  dr = '/xchip/cga/reference/annotation/db/ucsc/hg18';
elseif strcmpi(build,'hg19')
  dr = '';
elseif strcmpi(build,'hg19gencode')
  dr = '/xchip/cga/reference/annotation/db/ucsc/hg19';
elseif strcmpi(build,'mm9')
  dr = '/xchip/cga/reference/annotation/db/ucsc/mm9';
else
  error('unknown build %s',build);
end

fname = [dr 'cytoBand.txt'];
C = load_struct_noheader(fname,'%s%f%f%s%s');
if strcmp(C.col1{1},'chr')
  C = load_struct(fname,'%s%f%f%s%s');
else
  C = rename_fields(C,coln(1:5),{'chr','start','end','band','stain'});
end

C.chr = convert_chr(C.chr,build);
