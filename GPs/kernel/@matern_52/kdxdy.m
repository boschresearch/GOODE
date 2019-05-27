function dKd = kdxdy (k, x, y)
% KDXDY - Computes partial derivative covariance

% kdxdy.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-01-18
% Version: 0.1

x = x(:);
y = y(:);

  d = sqrt(5)/k.l * delta(k, x, y);
  
  D = - 25/(3*k.l^4) * exp(-d);
  
  dKd = NaN(size(x,1),size(y,1),size(x,2),size(y,2));
  for i=1:size(x,2)
    for j=1:size(y,2)
      dKd(:,:,i,j) = D .* bsxfun(@minus, x(:,i), y(:,i).') ...
                       .* bsxfun(@minus, x(:,j), y(:,j).');
                     
      if i==j
        dKd(:,:,i,j) = dKd(:,:,i,j) + 5/(3*k.l^2) * exp(-d) .* (1 + d);
      end
    end
  end
 
end