function len = load_chrlen(P)
% load_chrlen(method)
%
% method=='cytobands': based on cytobands definitions
%
% method=='chrfiles' (default) : based on filesize of hg18 chromosome files
%
% Mike Lawrence 2008-9


if ~exist('P','var'), P=[]; end

if ischar(P)
  build=P;
  P=[];
end

if ~exist('build','var')
  build = [];
end

P = impose_default_value(P,'build',build);

if isempty(P.build)
   fprintf('Assuming hg19\n');
   P.build = 'hg19';
%  P.build = 'chrfiles';
end

if strcmpi(P.build,'cytobands')

  fprintf('Assuming hg18\n');
  B=load_cytobands('hg18');
  len = zeros(24,1);
  for i=1:24
    idx = find(B.chr==i);
    len(i) = max(B.end(idx));
  end

elseif strcmpi(P.build,'chrfiles') ||strcmp(P.build,'hg18')

  if strcmpi(P.build,'chrfiles')
    fprintf('Assuming hg18\n');
  end
 
  len = [...    % filesize of hg18 chromosome files
     247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,...
     140273252,135374737,134452384,132349534,114142980,106368585,100338915,...
     88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954]';

%elseif strcmpi(P.build,'hg19')
%  len = [...
%       249250621,...
%243199373,...
%198022430,...
%191154276,...
%180915260,...
%171115067,...
%159138663,...
%146364022,...
%141213431,...
%135534747,...
%135006516,...
%133851895,...
%115169878,...
%107349540,...
%102531392,...
%90354753,...
%81195210,...
%78077248,...
%59128983,...
%63025520,...
%48129895,...
%51304566,...
%155270560,...
%59373566,...
%]';
%
%elseif strcmpi(P.build,'mm9')
else
  B=load_cytobands(P.build);
  ct = get_chrcount(P.build);
  len = zeros(ct,1);
  for i=1:ct
    idx = find(B.chr==i);
    len(i) = max(B.end(idx));
  end


%else
%  error('Unknown P.build');
end
