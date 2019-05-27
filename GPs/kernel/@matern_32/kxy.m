function K = kxy (k, x, y)
% KXY - Computes the Gram matrix between data x and y

  D = sqrt(3) * delta(k, x, y) / k.l;
  
  K = (1 + D) .* exp(- D);

end