% This source code is from GP_ODE_SOLVER
%  http://www.probabilistic-numerics.org/GP_ODE_Solver.zip
% Copyright (c) 2013-2014 Michael Schober, Philipp Hennig, Søren Hauberg
% This source code is licensed under the MIT license found in the
% 3rd-party-licenses.txt file in the root directory of this source tree.

function [ gp_kernel ] = se_kernel( ell, alpha )
% Construct a class of a squared exponential kernel defined as
%   k(x, y) = alpha * exp( -0.5 (x - y)^2 / ell^2 )
    
    % Prepare the struct.
    gp_kernel = struct();
    
    % Keep the lengthscale parameter.
    gp_kernel.ell = ell;
    gp_kernel.ell2 = ell^2;
    
    % Keep the kernel amplitude parameter.
    if( nargin > 1)
        gp_kernel.alpha = alpha;
    else
        gp_kernel.alpha = 1.0;
    end
    
    gp_kernel = class(gp_kernel, 'se_kernel'); 
end

