%  Purpose: Plot with relative error over Newton iterations and npoints
%  Author: John
%  Creation: 2019-01-20
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

startup;

%%
problem_index = 20;
bvp = str2func(['bvpT' int2str(problem_index)]);

% hyperparameters
maxit = 50;
abstol = 1e-5;
kernel = @se_kernel; %@matern_52; %
npoints_vec = [5,9,13,19,31,41,51,71,91,111];
N = length(npoints_vec);
ee = 1e-1;
M = 40;
%ls_fac_vec = linspace(1.5,15,M);
ls_fac_vec = logspace(log10(1.5),log10(15),M);

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-12);
solver = @bvp5c;

[probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvp(ee);
[problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();

% compute reference
sol_ref = solver(odefun,bcfun,solinit, options);

% iterate over npoints_vec
for kk = 1:N
    npoints = npoints_vec(kk);
    data{kk} = gpbvp_gridsearch(bvp, ee, npoints, ls_fac_vec, kernel, abstol, maxit);
    best_ls_fac(kk) = data{kk}.best_ls_ref_ind;
end

%% plots

%% error in each iteration for several npoints
themaxit = 7;
abstol = 1e-14;
error_it = NaN(themaxit,N);

for jj = 1:themaxit
    maxit = jj;
    for kk = 1:N
        npoints = npoints_vec(kk);
        ind = best_ls_fac(kk);
        ls_fac = ls_fac_vec(ind);
        [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'maxit', maxit, 'abstol', abstol);
        try % exact solution
            sol_ref_y = esolus(sol.x);
            temp = sol_ref_y(1,:);
        catch e % computed reference solution
            sol_ref_y = deval(sol_ref, sol.x);
            fprintf(2,'There was an error! The message was: %s\n',e.message);
        end 
        error_it(jj,kk) = norm(sol.y(:) - sol_ref_y(:))/norm(sol_ref_y(:));
    end
end


%% plot lines

close all
figure()
for kk = 1:N
    semilogy(1:themaxit,error_it(:,kk))
    hold on
end
legend
ylabel('relative error')
xlabel('iteration')

%% plot 3D surface

figure()
colormap winter
[X,Y] = meshgrid(npoints_vec,1:themaxit);
Z = error_it;
surf(X,Y,Z)
set(gca,'ZScale','log')
xlim([npoints_vec(1),npoints_vec(end)])
ylim([1,themaxit])
zlim([min(error_it(:)),max(error_it(:))])
xlabel('$N$', 'Interpreter', 'latex')
ylabel('Iteration')
zlabel('RelErr')
% colorbar
set(gca,'colorscale','log')
legend('hide')


