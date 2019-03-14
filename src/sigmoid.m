function y = sigmoid(x,inflection_point,sigmoidicity)
% y = sigmoid(x,inflection_point,sigmoidicity)

if ~exist('b','var'), b=1; end

x = double(x);
inflection_point = double(inflection_point);
sigmoidicity = double(sigmoidicity);

tmp = x;
tmp = tmp ./ inflection_point;
tmp = tmp .^ sigmoidicity;
tmp = tmp ./ (1+tmp);
y = tmp;




