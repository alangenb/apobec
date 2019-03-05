function z = matrix2lego(in)
% input is a 4-dimensional matrix of (left,base,right,newbase)
% output is an 8x12 matrix ready to serve as heights in the LEGO plot

if length(size(in))==4 && all(size(in)==[4 4 4 4])
  % left,base(ACGT),right,newbase
  % do strand-collapse
  for left=1:4, for right=1:4, for newbase=1:4
    in(left,[1 2],right,newbase) = in(left,[1 2],right,newbase) + in(5-right,[4 3],5-left,5-newbase);
  end,end,end
  in = in(:,1:2,:,:);
end


if ~(length(size(in))==4 && all(size(in)==[4 2 4 4]))
  error('unknown dimensions');
end

% left,base(AC),right,newbase

lego_classes = {'C->G','C->A','C->T','A->T','A->C','A->G'};
lego_bases = 'TCAG';
bases = 'ACGT';

z = zeros(8,12);
cx=1; cy=1;
for ci=1:6
  for left=1:4
    for right=1:4
      base_idx = lego_classes{ci}(1)==bases;
      newbase_idx = lego_classes{ci}(4)==bases;
      left_idx = lego_bases(left)==bases;
      right_idx = lego_bases(right)==bases;
      y = cy+left-1;
      x = cx+right-1;
      z(y,x,:) = in(left_idx,base_idx,right_idx,newbase_idx);
  end,end
  cx=cx+4;
  if cx>9, cx=1;cy=cy+4; end
end


