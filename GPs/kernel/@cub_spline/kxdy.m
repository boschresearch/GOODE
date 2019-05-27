function S = kxdy (k, x, y)
% KXDY - Computes the partial derivate of the covariance matrix with
% respect to the second input

% kxdy.m
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
  d = bsxfun(@minus, x, y.') .* d;  
  S = (1 + k.b) .* x + k.b .* ( - d - y.'.^2);
end