% This source code is from GP_ODE_SOLVER
%  http://www.probabilistic-numerics.org/GP_ODE_Solver.zip
% Copyright (c) 2013-2014 Michael Schober, Philipp Hennig, Søren Hauberg
% This source code is licensed under the MIT license found in the
% 3rd-party-licenses.txt file in the root directory of this source tree.

function [ Kddxddy] = kddxddy( gp_kernel, x, y )
% This evaluates the derivative of the kernel ddx and ddy.
    Kddxddy = (delta(gp_kernel, x, y).^4 - 6 * (delta(gp_kernel, x, y).^2) / gp_kernel.ell2 + 3 / (gp_kernel.ell2^2)) .* kxy(gp_kernel, x, y);
end