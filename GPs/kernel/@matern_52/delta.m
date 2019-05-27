function d = delta (k, x, y)
% DELTA - Computes the distance matrix of the kernel

% delta.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-05-19
% Version: 0.1

x = x(:);
y = y(:);
d = sqrt(bsxfun(@plus, dot(x,x,2), dot(y,y,2).') - 2*(x*y.'));

end
