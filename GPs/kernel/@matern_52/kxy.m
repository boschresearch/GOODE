function K = kxy (k, x, y)
% KXY - Computes the Gram matrix between data x and y

% kxy.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-05-19
% Version: 0.1

  D = sqrt(5) * delta(k, x, y) / k.l;
  
  K = (1 + D + D.^2 / 3) .* exp(- D);

end