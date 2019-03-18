function randinit(randseed)

if ~exist('randseed','var'), randseed=1234; end

rand('twister',randseed);
randn('seed',randseed);


