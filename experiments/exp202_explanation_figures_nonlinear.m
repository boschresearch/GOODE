%  Purpose: Create figures for explanation - nonlinear iterations
%  Author: John
%  Creation: 2019-01-12
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

%% Nonlinear solver

problem_index = 20;
bvp = str2func(['bvpT' int2str(problem_index)]);

% hyperparameters
kernel = @se_kernel ; %@matern_52; %
npoints = 31;
maxit = 50;
ee = 1e-1;
M = 20;
ls_fac_vec = linspace(.1,10,M);

error = NaN(1,M);
loglike = NaN(1,M);
sigma2 = NaN(1,M);

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-12);
solver = @bvp5c;

[probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvp(ee);
[problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();

% compute reference
sol_ref = solver(odefun,bcfun,solinit, options);

% iterate over ls_fac_vec
for jj = 1:M
    ls_fac = ls_fac_vec(jj);
    try
        disp(['ls_fac = ' num2str(ls_fac)])
        [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'MaxIt', maxit);
        disp(info)
        try % exact solution
            sol_ref_y = esolus(sol.x);
            temp = sol_ref_y(1,:);
        catch e % computed reference solution
            sol_ref_y = sol_ref.y;
            temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x);
            fprintf(2,'There was an error! The message was: %s\n',e.message);
        end
        error(jj) = norm(sol.y(1,:) - temp)/norm(temp);
        sigma2(jj) = mean(sol.sd(1,:));
        % check if loglike tends to -Inf
        loglike_max = max(loglike(1:jj-1));
        rel_change = abs(loglike_max-sol.loglike)/abs(loglike_max);
        if rel_change> 1e4
            disp(rel_change)
            loglike(jj) = -Inf;
            % break
        else
            loglike(jj) = sol.loglike;
        end
    catch e %e is an MException struct
        fprintf(2,'The identifier was: %s\n',e.identifier);
        fprintf(2,'There was an error! The message was: %s\n',e.message);
    end
end

%%

figure()
semilogy(ls_fac_vec, error)
xlabel('ls\_fac')
title('Lengthscale-factor vs. rel. error')
figure()
plot(ls_fac_vec, loglike)
xlabel('ls\_fac')
title('Lengthscale-factor vs. loglikelihood')

%% best solution according to logike

[~, ii] = max(loglike);
ls_fac = ls_fac_vec(ii);
disp(['Best ls_fac (log_like)= ' num2str(ls_fac)])

%% best solution according to reference
figure()
[~, ii] = min(error);
ls_fac = ls_fac_vec(ii);
disp(['Best ls_fac = ' num2str(ls_fac)])
[sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'MaxIt', maxit);
disp(info)
try
    sol_ref_y = esolus(sol.x);
    sol_ref_x = sol.x;
    temp = sol_ref_y(1,:);
catch
    sol_ref = bvp4c(odefun,bcfun,solinit);
    sol_ref_y = sol_ref.y;
    sol_ref_x = sol_ref.x;
    temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x)';
end

for ii=1:m
    subplot(1,m,ii)
    plot(sol_ref_x, sol_ref_y(ii,:),'r')
    hold on
    plot(sol.x, sol.y(ii,:), 'k-');
    hold on
    plot(sol.x, sol.y(ii,:)+3*sol.sd(ii,:), 'k:');
    hold on
    plot(sol.x, sol.y(ii,:)-3*sol.sd(ii,:), 'k:');
    xlabel('$t$')
    ylabel(['$y_' num2str(ii) '(t)$'])
    legend('exact', 'mean', '3sd')
end

figure()
for ii=1:m
    subplot(1,m,ii)
    plot(sol.x, sol.sd(ii,:), 'k');
    xlabel('$t$')
    ylabel(['$sd(y_' num2str(ii) '(t))$'])
end

%% local error, see Proposition 4.2 in [Cockayne 2016]
for ii=1:m
    figure()
    semilogy(sol_ref_x, abs(sol_ref_y(ii,:)-sol.y(ii,:)), 'Color', color_vec(end-1,:), 'Linestyle', ':', 'LineWidth', 1.0)
    hold on
    plot(sol.x, sol.sd(ii,:), 'Color', color_vec(end-1,:), 'Linestyle', '-.', 'LineWidth', 1.0);
    plot(sol.x, sol.sd(ii,:) * norm(sol_ref_y(ii,:)), 'Color', color_vec(end-1,:), 'Linestyle', '-', 'LineWidth', 1.0);
    xlabel('$t$', 'Interpreter','latex')
    ylabel(['$\sigma_' num2str(ii) '(t)$'], 'Interpreter','latex')
    legend('abs(exact-mean)', 'sd*norm(exact)')
    legend('hide')
    TikzFigure(['02_explanation_nonlin_gpbvp_sd_' num2str(ii) '.tex'])

end


%% pick single iterations and plot
maxit_vec = [1,2,3,4];
linestyle = {':', '-.', '--', '-', '-'};
locations = {'southwest', 'northwest'};

for ii=1:m
    fig = figure(20+ii);
    % set(fig,'defaultAxesColorOrder',[[1 1 1]; color_vec(end,:); color_vec(end,:)]);
    hold on
    plot(sol_ref_x, sol_ref_y(ii,:), 'LineWidth', 2.0, 'Color', color_vec(end-1,:), 'DisplayName', 'exact')
    xlabel('$t$', 'Interpreter','latex')
    ylabel(['$y_' num2str(ii) '(t)$'], 'Interpreter','latex')
    legend('Location', locations{ii})
end

for jj = 1:numel(maxit_vec)
    % compute solution for given MaxIt
    maxit = maxit_vec(jj);
    [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'MaxIt', maxit);
    
    for ii=1:m
        figure(20+ii)
        plot(sol.x, sol.y(ii,:), 'LineWidth', 1.0, 'DisplayName', ['It' num2str(maxit)], 'Color', color_vec(jj,:), 'linestyle', linestyle{jj});
        % plot sd - not very interesting changes a bit, but not much in
        % magnitude
%         figure(21)
%         hold on
%         subplot(1,m,ii)
%         plot(sol.x, sol.sd(ii,:), 'k:','DisplayName', ['It' num2str(maxit)]);
%         legend
    end
end

% standard deviation with second axis
% for ii=1:m
%     fig = figure(20+ii);
%     set(fig,'defaultAxesColorOrder',[[1 1 1]; color_vec(end,:); color_vec(end,:)]);
%     yyaxis right
%     plot(sol.x, sol.sd(ii,:), 'Color', color_vec(end,:), 'LineStyle', '-.');
%     %plot(sol_ref_x, abs(sol_ref_y(ii,:)-sol.y(ii,:)), 'Color', color_vec(end,:), 'LineStyle', ':')
%     ylabel(['$\sigma_' num2str(ii) '(t))$'], 'Interpreter','latex')
% end

for ii=1:m
    figure(20+ii)
    legend('hide')
    TikzFigure(['02_explanation_nonlin_gpbvp_' num2str(ii) '.tex'])
end



%% error in each iteration
themaxit = 10;
% npoints = 71;
error_it = NaN(themaxit);
for jj = 1:themaxit
    maxit = jj;
    [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'MaxIt', maxit, 'abstol', 1e-12);
    try % exact solution
        sol_ref_y = esolus(sol.x);
        temp = sol_ref_y(1,:);
    catch e % computed reference solution
        sol_ref_y = sol_ref.y;
        temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x);
        fprintf(2,'There was an error! The message was: %s\n',e.message);
    end
    
    error_it(jj) = norm(sol.y(1,:) - temp)/norm(temp);
    
end

figure()
semilogy(1:themaxit,error_it, 'k')
hold on
semilogy(sol.step, 'b')
semilogy(sol.residual_rel , 'g')
hold off
% semilogy(sol.rhs_norm)
legend('RelErr', 'step', 'RelRes')
ylabel('relative error')
xlabel('iteration')


