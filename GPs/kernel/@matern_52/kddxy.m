function ddK = kddxy (k, x, y)
% KDdXY - Computes mixed partial derivative covariance

% kddxy.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-12-02
% Version: 0.1

  d = sqrt(5)/k.l * delta(k, x, y);
  
  D = 25/(3*k.l^4) * exp(-d);
  
  ddK = NaN(size(x,1),size(y,1),size(x,2),size(x,2));
  for i=1:size(x,2)
    for j=1:size(x,2)
      ddK(:,:,i,j) = D .* bsxfun(@minus, x(:,i), y(:,i).') ...
                       .* bsxfun(@minus, x(:,j), y(:,j).');
      if i==j
        ddK(:,:,i,j) = ddK(:,:,i,j) - 5/(3*k.l^2) * exp(-d) .* (1 + d);
      end
    end
  end
 
end