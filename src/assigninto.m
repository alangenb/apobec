function A = assigninto(A,B,idx,flds1,flds2)
% same as mapinto, but with mapping already done
% TO DO: recode instead of just reformatting for mapinto

if ~exist('flds1','var'), flds1=fieldnames(B); end
if ~exist('flds2','var'), flds2=flds1; end

if isempty(A)
  dfname = 'dummy_field_34769348';
  A.(dfname) = nan(length(idx),1);
  flag = 1;
else
  flag = 0;
end

if slength(A)~=length(idx), error('A and idx should be same length'); end

f = 'idx_TEMP34567789';  % temporary fieldname
if isfield(A,f)||isfield(B,f), error(['please don''t use fieldname "' f '"']); end
B.(f) = as_column(1:slength(B));
A.(f) = nan(slength(A),1);
ii = find(~isnan(idx)&idx>=1&idx<=slength(B));
A.(f)(ii) = idx(ii);
A = mapinto(A,B,f,flds1,flds2);
A = rmfield(A,f);

if flag
  A = rmfield(A,dfname);
end



