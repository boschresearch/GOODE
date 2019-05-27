function dK = kdxy (k, x, y)
% KDXY - Computes mixed partial derivative covariance

x = x(:);
y = y(:);

  D = sqrt(3)/k.l * delta(k, x, y);
  
  dK = sqrt(3)/k.l * exp(-D) .* (- D);
end