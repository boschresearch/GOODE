function dK = getGradParams (k, x)
% GETGRADPARAMS - something

  D = sqrt(5) * delta(k, x, x);
  
  a = (1 + D ./ k.l + D.^2 / (3 * k.l^2));
  b = exp(-D ./ k.l);
  
  da = -D ./ k.l^2 - 2/3 * D.^2 ./ k.l^3;
  db = D ./ k.l^2 .* b;
  
  dK = reshape(da .* b + db .* a, [size(D), 1]);

end