function S = kxy (k, x, y)
% KXY - Computes the covariance matrix between inputs x and y

% kxy.m
% Author: David John
% Date: 2019-01-15
% -------------------------------------------------------------------
% Copyright (c) 2019 Robert Bosch GmbH
% All rights reserved.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
% 
% Authors: 
%    David John (david.john@de.bosch.com)
%    Michael Schober (michael.schober5@de.bosch.com)


  x = x(:);
  y = y(:);

  d = sqrt(bsxfun(@plus, dot(x,x,2), dot(y,y,2).') - 2*(x*y.'));
  d = d.^3 - bsxfun(@plus, x.^3, y'.^3); 
  S = 1 + (1 + k.b).* x*y' + k.b/3 .* d;
end