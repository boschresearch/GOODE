function dK = kdxy (k, x, y)
% KDXY - Computes mixed partial derivative covariance

% kdxy.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-05-19
% Version: 0.1

x = x(:);
y = y(:);

  D = sqrt(5)/k.l * delta(k, x, y);
  
  D = - 5/(3*k.l^2) * exp(-D) .* (1 + D);
  
  dK = NaN(size(x,1),size(y,1),size(x,2));
  for i=1:size(x,2)
    dK(:,:,i) = D .* bsxfun(@minus, x(:,i), y(:,i).');
  end
 
end