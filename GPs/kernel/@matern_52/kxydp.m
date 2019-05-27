function K = kxydp (k, x, y)
% KXYDP - Computes the derivative of the Gram matrix w.r.t. its parameters

% kxydp.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-02-24
% Version: 0.1

  D = sqrt(5) * delta(k, x, y) / k.l;
  
  K = 1/(3*k.l) * (D.^2 + D.^3) .* exp(- D);

end