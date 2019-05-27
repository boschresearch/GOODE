function [sol, info] = gpbvp_nonlinear_2(odefuns, bcfuns, dodefuns, dbcfuns, solinit, varargin)
%GPBVP  Solve two-point boundary value problems (BVP) using the Gaussian 
%       process formulation. The BVP is
%       y'(t) = f(t,y(t)), a <= t <= b, y \in R^D
%       g(y(a),y(b)) = 0.
%
% Inputs:
%   odefun  - a function handle that returns the ODE function f(t,y(t))
%   bcfun   - a function handle that returns the loss of the boundary cond.
%             g(y(a),y(b))
%   dodefun - a function handle that returns the Jacobian J_f(t,y(t))
%   dbcfuns - a function handle that returns the partial derivatives of
%             g(y(a),y(b)) with respect to both arguments
%   solinit - a structure containing an initial mesh solinit.x and an
%             initial guess solinit.y
%   varargin - can be used to pass additional options to the solver, see
%              below for more information
%                 alpha = 1.; % prior standard deviation
%                 ls_fac = 1.5; % factor for kernel lengthscale
%                 reg_fac = 1e-10; % regularization matrix for inversion
%                 npoints = 11; % number of points
%                 nout_fac = 10; % output points = npoints * noutfac
%                 kernel = @se_kernel; % default kernel
%                 MaxIt = 50; % max no of iterations
%                 AbsTol = 1e-6; % absolute stoping tolerance for Newton step
%
% Returns:
%   sol     - a structure containing the solution of the problem with sol.x
%             the prediction mesh and sol.y the approximate solution.
%   info    - some info
%
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
                                 

% TODO: check if derivatives are provided, else approximation is needed!!

% Set default values
alpha = 1.; % prior standard deviation
ls_fac = 1.5; % factor for kernel lengthscale
reg_fac = 1e-10; % regularization matrix for inversion
npoints = 11; % number of points
nout_fac = 10; % output points = npoints * noutfac
kernel = @se_kernel; % default kernel
MaxIt = 50; % max no of iterations
AbsTol = 1e-6; % absolute stoping tolerance for Newton step

if ~isempty(varargin)
    ii = 1;
    while ii <= length(varargin)
        if ischar(varargin{ii})
            switch lower(varargin{ii})
                case 'kernel'
                    kernel = varargin{ii+1};
                case 'ls_fac'
                    ls_fac = varargin{ii+1};
                case 'abstol'
                    AbsTol = varargin{ii+1};
                case 'maxit'
                    MaxIt = varargin{ii+1};
                case 'reg_fac'
                    reg_fac = varargin{ii+1};
                case 'alpha'
                    alpha = varargin{ii+1};
                case 'npoints'
                    npoints = varargin{ii+1};
                    % modify solinit
                    solinit = bvpinit(linspace(solinit.x(1),solinit.x(end),npoints),solinit.yinit);
                case 'nout_fac'
                    nout_fac = varargin{ii+1};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
        end
    ii = ii+1;    
    end
end

% initialize
sol = solinit;
[D, n] = size(solinit.y);
T = solinit.x;
DT = [T(1), T(end)];
N = nout_fac*n;
t = linspace(DT(1), DT(2), N); % prediction time points
% y0 = solinit.y(:,1);
% y1 = solinit.y(:,end);
% yab = -bcfuns(y0,y1);
[Ba, Bb] = dbcfuns();% e.g. Ba = [1, 0; 0, 0]; Bb = [0, 0; 1, 0];
muT_jj = reshape(solinit.y,D*n,1);

% ----- The kernel function
ell = median(diff(solinit.x))*ls_fac; % TODO heuristic for the length-scale
ker = kernel(ell);
V   = eye(D); % + ones(D)*0.5; % inter-dimensional cross-covariance

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT), kxy(ker, DT, T), kxdy(ker, DT, T);
     kxy(ker, T, DT), kxy(ker, T, T), kxdy(ker, T, T);
     kdxy(ker, T, DT), kdxy(ker, T, T), kdxdy(ker, T, T)];

F = [kxy(ker, t, DT), kxy(ker, t, T), kxdy(ker, t, T)];
FT = [kxy(ker, T, DT), kxy(ker, T, T), kxdy(ker, T, T)];

F2 = [kxy(ker, DT, t); kxy(ker, T, t); kdxy(ker, T, t)];
% F2T = [kxy(ker, DT, T); kxy(ker, T, T); kdxy(ker, T, T)];

Ktt = kxy(ker, t, t);
% KTT = kxy(ker, T, T);

% Iteration (Newton)
jj = 1;
step = 1e10;

while jj <= MaxIt && step >= AbsTol
    muT_jj_old = muT_jj;
    muT_jj_re = reshape(muT_jj,D,n);

    % calculate Y
    rhs_jj = odefuns(T, muT_jj_re);
    A = dodefuns(T, muT_jj_re); % Jacobian
    for ii = 1:n
        rhs_jj(:,ii) = rhs_jj(:,ii) - A(:,:,ii)*muT_jj_re(:,ii);
    end
    yab_jj = Ba*muT_jj_re(:,1) + Bb*muT_jj_re(:,end) - bcfuns(muT_jj_re(:,1),muT_jj_re(:,end));
    Y_jj = [yab_jj; rhs_jj(:)];
    
    % calculate H
    A_extended_jj = zeros(D*n); % if A const, then A_extended = kron(-1*eye(n),A)
    for ii = 1:n
        index = D*(ii-1)+1:D*ii;
        A_extended_jj(index,index) = A(:,:,ii);
    end
    H_jj = [Ba, Bb, zeros(D, D*n), zeros(D, D*n); 
        zeros(D*n, D), zeros(D*n, D), -1*A_extended_jj, kron(eye(n),eye(D))];
    L_jj = chol(  H_jj*kron(G,V)*H_jj'  + reg_fac*eye(D*n+D));
    % conditional mean
    temp = (L_jj \ (L_jj' \ Y_jj));
    muT_jj = (kron(FT,V) * H_jj') * temp;
    % covariance matrix
    % kT_post_jj = kron(KTT,V) - (kron(FT,V) * H_jj') * (L_jj \ (L_jj' \ (H_jj* (kron(F2T,V)))));

    % damping
    eta = 1.0;
    muT_jj = eta * muT_jj + (1-eta) * muT_jj_old;
    
    % L2 norm of Newton step
    step = norm(muT_jj_old(:) - muT_jj(:));
    sol.step(jj) = step;
    sol.step_rel(jj) = step/norm(muT_jj(:));

    % residual r_k = dy_k+1/dt - J_f(t,y_k) y_k+1 - b_k
    muT_jj_re = reshape(muT_jj,D,n); % new 
    rhs_norm = norm(rhs_jj(:));
    sol.rhs_norm(jj) = rhs_norm;
    for ii = 1:n
        rhs_jj(:,ii) = rhs_jj(:,ii) + A(:,:,ii)*muT_jj_re(:,ii);
    end
    residual = odefuns(T,muT_jj_re) - rhs_jj;
    residual = norm(residual(:));
    % relative residual
    sol.residual_rel(jj) = residual/rhs_norm;
    % forcing sequence (if relative residual <= eta, then Newton locally quadratic convergent)
    c = 1;
    sol.eta(jj) = min(0.5, c*rhs_norm);
    
    % increment
    jj = jj+1;
end

% compute and save last result
mu = (kron(F,V) * H_jj') * temp;
k_post = kron(Ktt,V) - (kron(F,V) * H_jj') * (L_jj \ (L_jj' \ (H_jj* (kron(F2,V)))));

sol.x = t;
sol.y = reshape(mu, D,N);
sol.sd = reshape(sqrt(diag(k_post)), D,N);
sol.k = k_post;

% log marginal likelihood
L_cov = L_jj;
log_det = 2* sum(log(diag(L_cov))); % \eqiv log(prod(diag(L_cov).^2));
LY = L_cov' \ Y_jj;
sol.loglike = -0.5*(LY' * LY) - 0.5* log_det - 0.5*(D*n+D)*log(2*pi);

info = {'points', n, 'iterations',jj-1, 'abserrlaststep', step};

end