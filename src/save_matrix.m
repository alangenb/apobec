function save_matrix(M,outfile)
% Mike Lawrence 2009-12-28

nr = size(M,1);
nc = size(M,2);

f = fopen(outfile,'wt');
for r=1:nr
  if nc>1, fprintf(f,'%d\t',M(r,1:end-1)); end
  fprintf(f,'%d\n',M(r,end));
end
fclose(f);
 
