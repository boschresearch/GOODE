% This source code is from GP_ODE_SOLVER
%  http://www.probabilistic-numerics.org/GP_ODE_Solver.zip
% Copyright (c) 2013-2014 Michael Schober, Philipp Hennig, Søren Hauberg
% This source code is licensed under the MIT license found in the
% 3rd-party-licenses.txt file in the root directory of this source tree.

function [ K ] = kxy( gp_kernel, x, y )
% This function evaluates the kernel.
    dist2 = bsxfun(@minus, x(:), y(:)').^2; % N x N
    K = gp_kernel.alpha * exp( -(0.5 / gp_kernel.ell2) * dist2 );

end

