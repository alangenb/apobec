function [ratio sd] = ratio_and_sd(n,N)

ratio = n./N;

sn = n .^ 0.5;
sN = N .^ 0.5;

sd = ratio .* ((sn./n).^2+(sN./N).^2).^0.5;


