function [val ttz] = get_from_fwb(varargin)
% get_from_fwb(fwb,chr,pos)
% get_from_fwb(fwb,chr,pos,build)
% get_from_fwb(fwb,chr,pos,build,fwi)    % where fwi is shared index for all members of fwb
% get_from_fwb(fwb,chr,start,end)
% get_from_fwb(chr,pos,fwb)
% get_from_fwb(chr,start,end,fwb)
% get_from_fwb(chr,start,end,fwb,build,fwi)
% etc.

% make sure jar is on classpath
javaclasspath('/xchip/cga/reference/mutsig_params/FixedWidthBinary.jar')

% FLEXIBLE ARGUMENT PARSE

args=[]; args.val = as_column(varargin);
args.can_be_fwb = false(nargin,1);
args.can_be_chr = false(nargin,1);
args.can_be_pos = false(nargin,1);
args.can_be_build = false(nargin,1);
args.can_be_fwi = false(nargin,1);
for i=1:slength(args)
  x = args.val{i};
  if ischar(x) && ~isempty(x)
    args.can_be_fwb(i) = true;
    args.can_be_fwi(i) = true;
    if strncmpi(x,'hg',2), args.can_be_build(i) = true; end
    if x(1)~='/', args.can_be_build(i) = true; end
    % possibly modify to check for if it's a build_dir
  elseif isnumeric(x)
    args.can_be_pos(i) = true;
    if ismember(x,[18 19 36 37])  % THIS MIGHT CAUSE PROBLEMS
      args.can_be_build(i) = true;
    end
    if ~any(x>50), args.can_be_chr(i) = true; end
  elseif iscellstr(x)
    if ~isempty(grep('^/',x))
      args.can_be_fwb(i) = true;
    else
      z  = str2double(x);
      if ~any(isnan(z)), args.can_be_pos(i) = true; end
      c = convert_chr(x);
      if mean(isnan(c))<0.1, args.can_be_chr(i) = true; end
    end
  end
end

% EXTRACT: chr pos   (OR)   chr start end
idx = find(args.can_be_chr,1);
if isempty(idx)
  error('requires chr argument');
else
  chr = args.val{idx};
  args = reorder_struct_exclude(args,idx);
end
idx = find(args.can_be_pos);
if length(idx)==1
  pos = args.val{idx(1)};
  args = reorder_struct_exclude(args,idx);
elseif length(idx)==2
  st = args.val{idx(1)};
  en = args.val{idx(2)};
  args = reorder_struct_exclude(args,idx);
else
  % (can't interpret)
end

% EXTRACT: fwb build fwi
idx = find(args.can_be_fwb,1);
if isempty(idx)
  error('requires fwb argument');
else
  fwb = args.val{idx};
  args = reorder_struct_exclude(args,idx);
end
idx = find(args.can_be_build,1);
if isempty(idx)
  % (ok)
else
  build = args.val{idx};
  args = reorder_struct_exclude(args,idx);
end
idx = find(args.can_be_fwi,1);
if isempty(idx)
  % (ok)
else
  fwi = args.val{idx};
  args = reorder_struct_exclude(args,idx);
end

% make sure no extraneous arguments
if slength(args)>=1
  error('unable to parse arguments');
end

% FUNCTION EXECUTION

if ischar(fwb), fwb={fwb}; end

if ~iscellstr(fwb), error('fwb should be a string filename or cellstr of filenames'); end

if exist('fwi','var') && iscellstr(fwi), error('multiple FWIs not supported'); end  % just loop over them

% make sure all FWBs+FWIs exist
%for i=1:length(fwb)
%  if length(fwb{i})<4 || ~strcmp(fwb{i}(end-3:end),'.fwb'), error('fwb should have "fwb" extension'); end
%end
demand_file(fwb);
if ~exist('fwi','var')
  demand_file(regexprep(fwb,'.fwb$','.fwi'));
else
  demand_file(fwi);
end

if ~exist('build','var')
  tmp = parse(fwb,'/db/([^/]+)/[^/]+/?$','build');
  if isempty(tmp.build)
    tmp.build = {''};  % leave it blank for convert_chr
  end
  build = tmp.build{1};
end
  
if ~isnumeric(chr), chr = convert_chr(chr,build); end

if exist('pos','var')
  if ~isnumeric(pos), pos = str2double(pos); end
  if length(chr)==1 && length(pos)>1
    chr = chr*ones(size(pos));
  end
  tot_query_len = length(pos);
else % st+en
  if ~isnumeric(st), st = str2double(st); end
  if ~isnumeric(en), en = str2double(en); end
  if length(st)~=length(en), error('length(st)~=length(en)'); end
  if length(chr)==1 && length(st)>1
    chr = chr*ones(size(st));
  end
  tot_query_len = sum(en-st+1);
end

if length(chr)>1e5
  fprintf('WARNING:  Long query, may take some time.\n');
  %   fprintf('WARNING:  Long query, may take some time.  Will print percentage indicator.\n');
  %  fprintf('WARNING:  get_from_fwb may not be the most efficient solution for very long position lists.\n');
  %  fprintf('          try reading from .mat files instead, if they exist.\n');
  %  fprintf('Continue?  (type "dbcont")\n');
  %  keyboard
end

if length(fwb)==1           % SINGLE-FILE MODE

  if exist('fwi','var')
    F = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwb{1},fwi);
  else
    F = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwb{1});
  end
  if exist('pos','var')
    if length(chr)>1e5
      %%%%%%%%%%%%%% BLOCK MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fprintf('Loading in block mode:');
      [cu ci cj] = unique(chr);
      val = nan(length(pos),1);
      for i=1:length(cu), fprintf(' %d/%d',i,length(cu));
        chr = cu(i);
        idx = find(cj==i & ~isnan(pos));
        posi = pos(idx);
        chrmin = min(posi);
        chrmax = max(posi);
        valrange = double(F.get(chr,chrmin,chrmax));
        val(idx) = valrange(posi-chrmin+1);
      end,fprintf('\n');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
      val = double(F.get(chr,pos));
    end
  else % st+en
    val = double(F.get(chr,st,en));
  end
  F.close();

  idx = find(val==-1);
  if ~isempty(idx)
    fprintf('%d/%d values in FWB were undefined (-1): these have been changed to NaN in the returned array.\n',...
            length(idx),tot_query_len);
    val(idx) = nan;
  end

else                        % MULTI-FILE MODE
 
  val = nan(tot_query_len,length(fwb));
  ttz = nan(length(fwb),1);

  if exist('fwi','var')
    F = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwb{1},fwi);
  else
    F = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwb{1});
  end
  step=1; if length(fwb)>30, step=10; end;  if length(fwb)>300, step=100; end
  fprintf('Loading from files: ');
  for i=1:length(fwb), if ~mod(i,step), fprintf('%d/%d ',i,length(fwb)); end
    tic

    if i>1
      if exist('fwi','var')    % using shared index
        F.switchFile(fwb{i});
      else
        F.close();
        F.open(fwb{i});
      end
    end

    if exist('pos','var')
      val(:,i) = double(F.get(chr,pos));
    else % st+en
      val(:,i) = double(F.get(chr,st,en));
    end
    ttz(i)=toc;

    idx = find(val(:,i)==-1);
    if ~isempty(idx)
      fprintf('%d/%d values in FWB were undefined (-1): these have been changed to NaN in the returned array.\n',...
              length(idx),tot_query_len);
      val(idx,i) = nan;
    end

  end
  fprintf('\n');
  F.close();

end

