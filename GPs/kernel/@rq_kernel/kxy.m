function S = kxy (k, x, y)
% KXY - Computes the covariance between inputs u and v of kernel k
%
% Inputs:
%   k - kernel
%   x - data matrix of first argument (row observations, columns data)
%   y - data matrix of second argument (row observations, columns data)

% kxy.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-04-28
% Version: 0.1

x = x(:);
y = y(:);

  % d = pdist2(u,v).^2; % see trick below - much faster!
  d = bsxfun(@plus, dot(x,x,2), dot(y,y,2).') - 2*(x*y.');
  d = d ./ k.lambda^2;
  
  % S = (1 + d ./ (2*k.alpha)).^(-k.alpha); 
  S = exp(-k.alpha * log(1 + d ./ (2*k.alpha))); % twice as fast :-]

end