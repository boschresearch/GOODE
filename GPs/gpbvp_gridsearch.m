function [data] = gpbvp_gridsearch(bvp, ee, npoints, ls_fac_vec, kernel, abstol, maxit)
%gpbvp_gridsearch: helper function to perfom gridsearch for kernel parameter optimization
%
% Inputs:
%   bvp     - a function handle that returns the problem setting
%   ee      - BVP parameter
%   and more ...

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


[probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvp(ee);
[problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();

esol = 0;
M = numel(ls_fac_vec); 

error = NaN(1,M);
loglike = NaN(1,M);
sigma = NaN(1,M);

% compute reference
options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-12);
solver = @bvp5c;
sol_ref = solver(odefun,bcfun,solinit, options);

% iterate over ls_fac_vec
for jj = 1:M
    ls_fac = ls_fac_vec(jj);
    try
        %disp(['ls_fac = ' num2str(ls_fac)])
        [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'abstol', abstol, 'maxit', maxit);
        %disp(info)
        try % exact solution
            sol_ref_y = esolus(sol.x);
            esol = 1;
        catch e % computed reference solution
            sol_ref_y = deval(sol_ref, sol.x);
            fprintf(2,'There was an error! The message was: %s\n',e.message);
        end
        error(jj) = norm(sol.y(:) - sol_ref_y(:))/norm(sol_ref_y(:));
        sigma(jj) = mean(sol.sd(:));
        % check if loglike tends to -Inf
        loglike_max = max(loglike(1:jj-1));
        rel_change = abs(loglike_max-sol.loglike)/abs(loglike_max);
        if 0 %rel_change> 1e8
            disp(rel_change)
            loglike(jj) = -Inf;
        else
            loglike(jj) = sol.loglike;
        end
    catch e %e is an MException struct
        fprintf(2,'The identifier was: %s\n',e.identifier);
        fprintf(2,'There was an error! The message was: %s\n',e.message);
    end
end

data.error = error;
data.loglike = loglike;
data.sigma = sigma;
data.linear = Linear;
data.esol = esol;


% best solution according to logike
[~, ii] = max(loglike);
data.best_ls_loglike_ind = ii;
data.best_ls_loglike = ls_fac_vec(ii);

%% best solution according to reference
[~, ii] = min(error);
data.best_ls_ref_ind = ii;
data.best_ls_ref = ls_fac_vec(ii);

end