function dS = kdxdy (k, x, y)
% KDXDY - derivative of the kernel w.r.t. to both arguments
%
% Depending on the input space of the parameters this returns a
% multi-dimensional matrix, i.e., kdxy : R^D x R^D ->  R^D. Dimensions are
% returned as third matrix dimension, i.e., for matrix-valued inputs
% X in R^N1xD, y in % R^N2xD, kdxy returns a N1xN2xD tensor.

% kdxdy.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-04-28
% Version: 0.1

x = x(:);
y = y(:);

  D = size(x,2);

  dsq = bsxfun(@plus, dot(x,x,2), dot(y,y,2).') - 2*(x*y.');
  
  dsq2 = (1 + dsq / (2 * k.alpha * k.lambda^2)); %.^(-k.alpha-1);
  
  dS = NaN(size(x,1),size(y,1),D);
  for i=1:D
    dS(:,:,i) = dsq2.^(-k.alpha-2) .* (-k.alpha-1) .* dsq./(k.alpha*k.lambda^4) + dsq2.^(-k.alpha-1) ./(k.lambda^2);
  end
end