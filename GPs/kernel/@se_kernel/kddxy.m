% This source code is from GP_ODE_SOLVER
%  http://www.probabilistic-numerics.org/GP_ODE_Solver.zip
% Copyright (c) 2013-2014 Michael Schober, Philipp Hennig, Søren Hauberg
% This source code is licensed under the MIT license found in the
% 3rd-party-licenses.txt file in the root directory of this source tree.

function [ Kddxy ] = kddxy( gp_kernel, x, y )
% This evaluates the derivative of the kernel ddx.
    Kddxy = -kdxdy(gp_kernel, x, y);
end

