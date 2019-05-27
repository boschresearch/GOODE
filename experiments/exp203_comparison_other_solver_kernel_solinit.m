%  Purpose: Comparison of GP-BVP solver to standard solvers for all
%  problems, with SE kernel and mesh from bvp5c. Also plots for Paper
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

%%
start = 1;
n_problems = 33;
bvp_list = {start:n_problems};
for ii = start:n_problems
    bvp_list{ii} = str2func(['bvpT' int2str(ii)]);
end
nbvp = numel(bvp_list);

% hyperparameters
% npoints = 71; We use the points from bvp5c 
ee = 1e-1;
M = 40;
% ls_fac_vec = linspace(1.5,15,M);
ls_fac_vec = logspace(log10(1.5),log10(15),M);
kernel = @se_kernel;

% ls_fac_vec = logspace(log10(1.5),log10(150),M);
% kernel = @matern_52;

% npoints = 71;
% ls_fac_vec = (logspace(log10(1),log10(1.5),M)-1).^2; % uiuiui need to change ls_fac to param b
% kernel = @cub_spline;

% tolerance for reference solution if no exact solution provided
options_ref = bvpset('RelTol', 1e-9, 'AbsTol', 1e-12);
solver_ref = @bvp5c;
% tol for other solvers
options = bvpset('RelTol', 1e-5, 'AbsTol', 1e-8);
solvers = {@bvp4c, @bvp5c @tom, @bvptwp};
n_solvers = numel(solvers);


error = NaN(nbvp,n_solvers);
points = NaN(nbvp,n_solvers);
deltat = NaN(nbvp,n_solvers);

error_gp = NaN(nbvp,M);
loglike = NaN(nbvp,M);
deltat_gp = NaN(nbvp,1);
esolus_mask = zeros(nbvp,1);


for ii = start:nbvp
    disp(['problem ' num2str(ii)])
    [probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvp_list{ii}(ee);
    [problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();
    
    % compute reference
    sol_ref = solver_ref(odefun,bcfun,solinit, options_ref);
    % compute bvp5c basic and use mesh for gpbvp
    sol_ref_b = solver_ref(odefun,bcfun,solinit, options);
    solinit_b = bvpinit(sol_ref_b.x,solinit.yinit);

    % deval(sol_ref,x); % evaluates solution at desired points
    npoints_ref = size(sol_ref.x,2);
    disp(npoints_ref)
    % deltat_ref = min(diff(sol_ref.x));
    % solinit_ref = bvpinit(sol_ref.x,solinit.yinit);
    
    % iterate over ls_fac_vec
    for jj = 1:M
        ls_fac = ls_fac_vec(jj);
        try
            %disp(['ls_fac = ' num2str(ls_fac)])
            [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit_b, 'ls_fac', ls_fac, 'kernel', kernel);
            %disp(info)
            try % exact solution
                esol = esolus(sol.x);
                temp = esol(1,:);
                esolus_mask(ii) = 1;
            catch % computed reference solution
                temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x);
            end
            error_gp(ii,jj) = norm(sol.y(1,:) - temp)/norm(temp);
            
            % check if loglike tends to -Inf
            loglike_max = max(loglike(ii,1:jj-1));
            rel_change = (loglike_max-sol.loglike)/abs(loglike_max);
            if rel_change> 1e3
                disp(rel_change)
                loglike(ii,jj) = -Inf;
            else
                loglike(ii,jj) = sol.loglike;
            end
        catch e %e is an MException struct
            fprintf(2,'The identifier was: %s\n',e.identifier);
            fprintf(2,'There was an error! The message was: %s\n',e.message);
        end
    end
    deltat_gp(ii) = min(diff(sol.x));

    
    % other solvers
    for jj = 1:n_solvers
        try
            sol = solvers{jj}(odefun,bcfun,solinit,options);
            % deval(sol_ref,x); % evaluates solution at points x
            try 
                esol = esolus(sol.x);
                temp = esol(1,:);
            catch
                sol_ref_y = sol_ref.y;
                temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x);
            end
            error(ii,jj) = norm(sol.y(1,:) - temp)/norm(temp);
            points(ii,jj) = max(size(sol.x));
            deltat(ii,jj) = min(diff(sol.x));
        catch e %e is an MException struct
            fprintf(2,'The identifier was: %s\n',e.identifier);
            fprintf(2,'There was an error! The message was: %s\n',e.message);
        end
    end
    
end
%% Creat some nice plots 

startup;
marker = {'x', '+', '^', 's', 'o'};
our_color = 7;

% Select best according to loglike
[~, ind] = max(loglike');
linear_ind = sub2ind(size(loglike),1:n_problems, ind);
error_opt_loglike = error_gp(linear_ind);

% figure(6)
% kkk = 40
% plot(ls_fac_vec(1:kkk),loglike(23,1:kkk)')
% figure(7)
% semilogy(ls_fac_vec(1:kkk),error_gp(23,1:kkk)')

% Select best according to error
[~, ind] = min(error_gp');
linear_ind = sub2ind(size(error_gp),1:n_problems, ind);
error_gp_min = error_gp(linear_ind);

figure()
semilogy(error_opt_loglike, ':^', 'LineWidth', 1.5, 'Color', color_vec(our_color-1,:), 'DisplayName', 'GOODE loglike')
hold on
semilogy(error_gp_min, ':o', 'LineWidth', 1.5, 'Color', color_vec(our_color,:), 'DisplayName', 'GOODE')
ylabel('relative error')
xlabel('problem no.')
% legend('Location','best')
legend off
xticks(1:2:33)
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';
% title('RelErr of loglike optimum')

TikzFigure('03_comparison_loglike_vs_global_opt_solinit.tex')

%% plot comparison of all solver - relative error

figure()
for jj = 1:n_solvers
    semilogy(1:33,error(:,jj), marker{jj}, 'Color', color_vec(jj,:), 'LineWidth', 1., 'DisplayName', func2str(solvers{jj}))
    hold on 
end
semilogy(1:33,error_gp_min, 'o', 'Color', color_vec(our_color,:), 'LineWidth', 2., 'DisplayName', 'GOODE')
ylabel('relative error')
xlabel('problem no.')
%legend('Location','best')
legend off
xticks(1:2:33)
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';
%title('RelErr')

TikzFigure('04_comparison_all_solver_relerr_solinit.tex')


%% plot no points

figure()
for jj = 1:n_solvers
    semilogy(1:33,points(:,jj), marker{jj}, 'Color', color_vec(jj,:), 'LineWidth', 1., 'DisplayName', func2str(solvers{jj}))
    hold on 
end
semilogy(1:33,points(:,2), 'o', 'Color', color_vec(our_color,:), 'LineWidth', 2., 'DisplayName', 'GOODE')
ylabel('N', 'Interpreter','latex')
xlabel('problem no.')
%legend('Location','best')
legend off
xticks(1:2:33)
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';
%title('No. points')

TikzFigure('04_comparison_all_solver_points_solinit.tex')

% make table
mean(points)
std(points)


%% plot deltat

figure()
for jj = 1:n_solvers
    semilogy(1:33,deltat(:,jj), marker{jj}, 'Color', color_vec(jj,:), 'LineWidth', 1., 'DisplayName', func2str(solvers{jj}))
    hold on 
end
semilogy(1:33,deltat_gp, 'o', 'Color', color_vec(our_color,:), 'LineWidth', 2., 'DisplayName', 'GOODE')
ylabel('$h$', 'Interpreter','latex')
xlabel('problem no.')
% legend('Location','best')
% title('deltat')
legend off
xticks(1:2:33)
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';

TikzFigure('04_comparison_all_solver_deltat_solinit.tex')


%% GP solver error and loglike vs. ls_fac
%figure()
% HOW to plot a nice picture???
% imshow(log(error))

figure()
semilogy(ls_fac_vec,error_gp(1:18,:)')
xlabel('lsfac')
ylabel('relative error')
title('Linear')
legend()

figure()
semilogy(ls_fac_vec,error_gp(19:33,:)')
xlabel('lsfac')
ylabel('relative error')
title('Nonlinear')
legend()

figure()
plot(ls_fac_vec,loglike(1:18,:)')
xlabel('lsfac')
ylabel('loglike')
title('Linear')
legend()

figure()
plot(ls_fac_vec,loglike(19:33,:)')
xlabel('lsfac')
ylabel('loglike')
title('Nonlinear')
legend()
%% look at specific problems
problem_idx = 23;

figure()
semilogy(ls_fac_vec,error_gp(problem_idx,:))
xlabel('lsfac')

figure()
semilogy(ls_fac_vec,-loglike(problem_idx,:))
xlabel('lsfac')

