%  Purpose: Test GP-BVP solver for all problems from BVPTWP
%  Author: John
%  Creation: 2018-08-31
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

clear all;
close all;

%%
n_problems = 33;
bvp_list = {1:n_problems};
for ii = 1:n_problems
    bvp_list{ii} = str2func(['bvpT' int2str(ii)]);
end
nbvp = numel(bvp_list);

% parameters
npoints = 31;
ls_fac = 3;
ee = 1e-1;

error = zeros(nbvp,1);
linear_mask = zeros(nbvp,1);
esolus_mask = zeros(nbvp,1);
dimension = zeros(nbvp,1);

for ii = 1:nbvp
%     if ii == 24
%         continue
%     end
    disp(['problem ' num2str(ii)])
    [probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvp_list{ii}(ee);
    [problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();
    disp(['Linear: ' Linear])
    if strcmp(Linear,'on')
        linear_mask(ii) = 1;
    end
    try
        esolus(0.1);
        disp('esolus: on')
        esolus_mask(ii) = 1;
    catch
        disp('esolus: off')
    end
    dimension(ii) = m;
%     sol_ref = bvp4c(odefun,bcfun,solinit);
%     
%     [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'npoints', npoints, 'ls_fac', ls_fac);
%     disp(info)
%     temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x)';
%     error(ii) = norm(sol.y(1,:) - temp)/norm(temp);
end

figure()
plot(1:n_problems, linear_mask, 'b*')
hold on
plot(1:n_problems, esolus_mask, 'ro')
plot(1:n_problems, dimension, 'k^')
legend('Linear', 'Exact solution', 'dimension')
xlabel('problem no.')


% semilogy(error)
% ylabel('relative error')
% xlabel('problem no.')

