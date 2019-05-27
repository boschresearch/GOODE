function Kd = kxdy (k, x, y)
% KXDY - Computes mixed partial derivative covariance

% kxdy.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-01-18
% Version: 0.1

x = x(:);
y = y(:);

  D = sqrt(5)/k.l * delta(k, x, y);
  
  D = 5/(3*k.l^2) * exp(-D) .* (1 + D);
  
  Kd = NaN(size(x,1),size(y,1),size(y,2));
  for i=1:size(y,2)
    Kd(:,:,i) = D .* bsxfun(@minus, x(:,i), y(:,i).');
  end
 
end