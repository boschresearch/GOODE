function dKd = kdxdy (k, x, y)
% KDXDY - Computes partial derivative covariance

x = x(:);
y = y(:);

  D = sqrt(3)/k.l * delta(k, x, y);
  
  dKd = 3/(k.l^2) *(1-D).* exp(-D);
end