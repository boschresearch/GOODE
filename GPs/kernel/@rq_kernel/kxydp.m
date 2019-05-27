function K = kxydp (k, x, y)
% KXYDP - Computes the derivative of the Gram matrix w.r.t. its parameters
%
% Inputs:
%   k - kernel
%   x - data matrix of first argument (row observations, columns data)
%   y - data matrix of second argument (row observations, columns data)

% kxydp.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-02-25
% Version: 0.1

  % dsq = r^2 / l^2
  dsq = bsxfun(@plus, dot(x,x,2), dot(y,y,2).') - 2*(x*y.');
  dsq = dsq ./ k.lambda^2;
  
  % P = (1 + r^2/2*a*l^2)
  P = 1 + dsq ./ (2*k.alpha);
  
  K = NaN([size(dsq), 2]);
  % Kda
  K(:,:,1) = P.^(-k.alpha) .* (1/(2*k.alpha) * dsq ./ P - log(P));
  % Kdl
  K(:,:,2) = 1/k.lambda * dsq .* P.^(-k.alpha-1);

end