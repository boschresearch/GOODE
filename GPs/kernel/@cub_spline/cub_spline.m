function ker = cub_spline (b, l)
% CUB_SPLINE - Constructor for the cubic spline kernel object

% cub_spline.m
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


  ker = struct ();
  
  if nargin < 2
    l = 1.; % renormalization such that inputs are in [0,1]
  end
  
  ker.b = b; % step weights
  ker.l = l; 
  
  ker = class(ker, 'cub_spline');

end