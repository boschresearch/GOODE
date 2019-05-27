function ker = matern_32 (l)
% MATERN_32 - Constructor for a Matern 3/2 kernel object defined as
%   k(x, y) = (1 + D)*exp(-D)
%   D = sqrt(3) * delta(k, x, y) / k.l;

  ker = struct ();
  
  ker.l = l; % characteristic length scale
  
  ker = class(ker, 'matern_32');

end