function tidx = map_mutations_to_targets_fast_wrapper(M,T,quiet_flag)
% returns NaN for unmapped mutations

if ~exist('quiet_flag','var'), quiet_flag = false; end

if ~isfield(M,'chr') && isfield(M,'chr1'), M.chr=M.chr1; end
if ~isfield(M,'pos') && isfield(M,'pos1'), M.pos=M.pos1; end
if ~isfield(M,'pos') && isfield(M,'st'), M.pos=M.st; end

M = keep_fields(M,{'chr','pos'});
if ~isfield(T,'start') && ~isfield(T,'st') && isfield(T,'gene_start'), T.st=T.gene_start; end
if ~isfield(T,'end') && ~isfield(T,'en') && isfield(T,'gene_end'), T.en=T.gene_end; end
T = require_fields_with_convert(T,{'chr','st','en'},{'chr','start','end'});
M.chr = double(convert_chr(M.chr));
M.pos = double(M.pos);
T.chr = double(convert_chr(T.chr));
T.st = double(T.st);
T.en = double(T.en);

M.idx = (1:slength(M))';
T.idx = (1:slength(T))';

if ~quiet_flag, fprintf('[sort mutations]  '); end
M = sort_struct(M,{'chr','pos'});

if ~quiet_flag,fprintf('[sort targets]  '); end
T = sort_struct(T,{'chr','st'});

if ~quiet_flag,fprintf('[map]  '); end
M.tidx = map_mutations_to_targets_fast([M.chr M.pos],[T.chr T.st T.en]);
if ~quiet_flag,fprintf('\n'); end

tidx = nan(slength(M),1);

T.idx(end+1) = nan; M.tidx(M.tidx==-1) = length(T.idx);   % to handle unmapped mutations

tidx(M.idx) = T.idx(M.tidx);

