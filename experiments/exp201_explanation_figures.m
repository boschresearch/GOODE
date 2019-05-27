%  Purpose: Create figures for explanation - linear 
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

%% Linear solver

seed = 1;
rng(seed);

SAVE_FIGURES = true;

problem_index = 1;
bvp = str2func(['bvpT' int2str(problem_index)]);

% hyperparameters
kernel = @se_kernel ; %@matern_52; %
npoints = 5;
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
        [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'MaxIt', 1.0);
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
% [~, ii] = min(error);
ls_fac = 1.3; % ls_fac_vec(ii);
disp(['Best ls_fac = ' num2str(ls_fac)])
[sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel, 'MaxIt', 1.0);
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
    ylabel(['$y_' num2str(ii-1) '(t)$'])
    legend('exact', 'mean', '3sd')
end

figure()
for ii=1:m
    subplot(1,m,ii)
    plot(sol.x, sol.sd(ii,:), 'k');
    xlabel('$t$')
    ylabel(['$sd(y_' num2str(ii-1) '(t))$'])
end

%% GP regression

% modify solinit
solinit = bvpinit(linspace(solinit.x(1),solinit.x(end),npoints),solinit.yinit);

% data points
T = solinit.x(2:end-1);
DT = [solinit.x(1), solinit.x(end)];
yd = esolus(solinit.x);
sigma = 0.25 * median(diff(solinit.x)); % too much?
% prediction
t = sol.x;
ell = 0.75 * median(diff(solinit.x))*ls_fac;
[D, n] = size(solinit.y);
ker = kernel(ell);

theta = 0.01;

G = [kxy(ker, DT, DT), kxy(ker, DT, T);
     kxy(ker, T, DT), kxy(ker, T, T)];
G = theta * G;   
   
F = [kxy(ker, t, DT), kxy(ker, t, T);
     kdxy(ker, t, DT), kdxy(ker, t, T)];
F = theta * F;
   
F2 = F';
Ktt = [kxy(ker, t, t), kxdy(ker, t, t);
       kdxy(ker, t, t), kdxdy(ker, t, t)];
Ktt = theta * Ktt;     
     
L = chol(G + blkdiag(zeros(numel(DT)), sigma^2*eye(numel(T))));
Y_clean = [yd(1,1); yd(1,end); yd(1,2:end-1)'];
Y = Y_clean + [zeros(2,1); sigma * randn(numel(T),1)];

%{
G = kxy(ker, T, T);
F = kxy(ker, t, T);
F2 = kxy(ker, T, t);
Ktt = kxy(ker, t, t);
V   = eye(D);
L = chol(  kron(G,V)  + 1e-10*eye(D*n));
Y = reshape(yd,D*n,1);
%}
mu = F * (L \ (L' \ Y));
% mu = reshape(mu,D,numel(t));
mu = reshape(mu,numel(t),2)';
k_post = Ktt - F * (L \ (L' \ F2));
% sd = reshape(sqrt(diag(k_post)),D,numel(t));
sd = reshape(sqrt(diag(k_post)),numel(t),D)';

figure()
for ii=1:m
    subplot(1,m,ii)
    plot(sol_ref_x, sol_ref_y(ii,:),'r')
    hold on
    plot([DT(1), T, DT(end)], yd(ii,:),'r*')
    hold on
    plot(sol.x, mu(ii,:), 'k-');
    hold on
    plot(sol.x, mu(ii,:)+3*sd(ii,:), 'k:');
    hold on
    plot(sol.x, mu(ii,:)-3*sd(ii,:), 'k:');
    xlabel('t')
    ylabel(['y_' num2str(ii-1) '(t)'])
    legend('exact', 'data', 'mean', '3sd')
    legend('Location', 'best')
end

figure()
for ii=1:m
    subplot(1,m,ii)
    plot(sol.x, sd(ii,:), 'k');
    xlabel('$t$')
    ylabel(['$sd(y_' num2str(ii-1) '(t))$'])
end

%% Final plot (all in one with offset)

figure()

% gp_bvp solution
for ii=1:m
    subplot(1,m,ii)
    hold on
    plot(sol_ref_x, sol_ref_y(ii,:),'r')
    plot(sol.x, sol.y(ii,:), 'k-');
    plot(sol.x, sol.y(ii,:)+3*sol.sd(ii,:), 'k:');
    plot(sol.x, sol.y(ii,:)-3*sol.sd(ii,:), 'k:');
    xlabel('$t$')
    ylabel(['$y_' num2str(ii-1) '(t)$'])
    legend('exact', 'mean', '3sd')
end

% standard regression

offset = 0.5;
for ii=1:m
    subplot(1,m,ii)
    hold on
    plot(sol_ref_x, sol_ref_y(ii,:) + offset,'r')
    plot([DT(1), T, DT(end)], yd(ii,:) + offset,'r*')
    plot(sol.x, mu(ii,:) + offset, 'k-');
    plot(sol.x, mu(ii,:)+3*sd(ii,:) + offset, 'k:');
    plot(sol.x, mu(ii,:)-3*sd(ii,:) + offset, 'k:');
    xlabel('$t$')
    ylabel(['$y_' num2str(ii-1) '(t)$'])
    legend('exact', 'data', 'mean', '3sd')
    legend('Location', 'best')
end

%% Final plot (seperate plots)
locations = {'northeast', 'southeast'};
fac_sd = 10;
our_color = 7;
reg_color = 3;

% gp_bvp solution
ylimits = NaN(m,2);

for ii=1:m
    figure()
    hold on
    plot(sol_ref_x, sol_ref_y(ii,:),'LineWidth', 2.0, 'Color', color_vec(end-1,:))
    plot(sol.x, sol.y(ii,:),'LineWidth', 1.0, 'Color', color_vec(our_color,:), 'linestyle', '--' );
    plot(sol.x, sol.y(ii,:)+fac_sd*sol.sd(ii,:),'LineWidth', 1.0, 'Color', color_vec(our_color,:), 'linestyle', ':' );
    plot(sol.x, sol.y(ii,:)-fac_sd*sol.sd(ii,:),'LineWidth', 1.0, 'Color', color_vec(our_color,:), 'linestyle', ':' );
    xlabel('$t$', 'Interpreter','latex')
    xticks(solinit.x)
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'off';
    ylabel(['$y_' num2str(ii) '(t)$'], 'Interpreter','latex')
    legend('exact', 'mean', [num2str(fac_sd) '*sd'])
    if ii ==1
        scatter([DT(1), DT(end)], [yd(ii,1), yd(ii,end)], 200, 'ks')
        legend('exact', 'mean', [num2str(fac_sd) '*sd'], [num2str(fac_sd) '*sd'], 'bc')
    end
    legend('Location', locations{ii})
    legend('hide')
    ylimits(ii,:) = ylim();
    if SAVE_FIGURES
    TikzFigure(['01_explanation_lin_gpbvp_' num2str(ii) '.tex'])
    end
end

% standard regression
offset = 0.0;
for ii=1:m
    figure()
    hold on
    plot(sol_ref_x, sol_ref_y(ii,:) + offset, 'LineWidth', 2.0, 'Color', color_vec(end-1,:))
    % scatter(T, yd(ii,:) + offset, 200, 'ks' )
    if ii ==1
        scatter([DT, T], Y, 200, 'ks')
        legend('exact', 'mean', [num2str(fac_sd) '*sd'], [num2str(fac_sd) '*sd'], 'bc')
    end
    plot(sol.x, mu(ii,:) + offset, 'LineWidth', 1.0, 'Color', color_vec(reg_color,:), 'linestyle', '--' );
    plot(sol.x, mu(ii,:)+fac_sd*sd(ii,:) + offset, 'LineWidth', 1.0, 'Color', color_vec(reg_color,:), 'linestyle', ':' );
    plot(sol.x, mu(ii,:)-fac_sd*sd(ii,:) + offset, 'LineWidth', 1.0, 'Color', color_vec(reg_color,:), 'linestyle', ':' );
    xlabel('$t$', 'Interpreter','latex')
    xticks(solinit.x)
    ylabel(['$y_' num2str(ii) '(t)$'], 'Interpreter','latex')
    ylim(ylimits(ii,:));
    legend('exact', 'data', 'mean', [num2str(fac_sd) '*sd'])
    legend('Location', locations{ii})
    legend('hide')
    if SAVE_FIGURES
    TikzFigure(['01_explanation_lin_gpreg_' num2str(ii) '.tex'])
    end
end