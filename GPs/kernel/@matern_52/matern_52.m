function ker = matern_52 (l)
% MATERN_52 - Constructor for a Matern 5/2 kernel object defined as
%   k(x, y) = (1 + D + 1/3*D^2)*exp(-D)
%   D = sqrt(5) * delta(k, x, y) / k.l;

% matern_52.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-05-19
% Version: 0.1

  ker = struct ();
  
  ker.l = l; % characteristic length scale
  
  ker = class(ker, 'matern_52');

end