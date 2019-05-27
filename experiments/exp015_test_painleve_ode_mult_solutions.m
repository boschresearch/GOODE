%  Purpose: Check problem 34 - Painleve ODE with 2 two solutions.
%  Need to change initialisation to find both solutions
%  init solution 1 is 0
%  init solution 2 is linear between (x,y)=(0,-3) and (10,3) for y_1 and 0
%  for y_2

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
ee = 1e-0;
[probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvpT34(ee);
[problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();

%% plot settings
startup;
our_color = 7;
bvp5c_color = 2;

%%
maxit = 50;
abstol = 1e-4;
M = 40;
npoints = 31;
kernel = @se_kernel;  % @matern_52; %
ls_fac_vec = logspace(0,log10(15),M);
error = NaN(2,M);
loglike = NaN(2,M);

% init for solution 1
init{1} = [0,0; 0,0];
% init for solution 2
init{2} = [-3,3; 0,0];

options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-6); % default 'RelTol', 1e-3, 'AbsTol', 1e-6


for kk =1:2 % iterate over both solutions
    % modify solinit
    t = linspace(solinit.x(1),solinit.x(end),npoints);
    yinit = zeros(m,npoints);
    for ii =1:m
        yinit(ii,:) = interp1([solinit.x(1),solinit.x(end)], init{kk}(ii,:), t);
    end
    solinit.x = t;
    solinit.y = yinit;
    
    % reference solution
    sol_ref = bvp5c(odefun,bcfun,solinit, options);
    disp(sol_ref.stats)
    
    for ii = 1:M
        ls_fac = ls_fac_vec(ii);
        try
            disp(['ls_fac = ' num2str(ls_fac)])
            [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit,  'ls_fac', ls_fac, 'kernel', kernel, 'maxit', maxit, 'abstol', abstol);
            disp(info)
            sol_ref_y = deval(sol_ref, sol.x);
            error(kk,ii) = norm(sol.y(:) - sol_ref_y(:))/norm(sol_ref_y(:));
            loglike(kk,ii) = sol.loglike;
        catch e %e is an MException struct
            fprintf(2,'The identifier was: %s\n',e.identifier);
            fprintf(2,'There was an error! The message was: %s\n',e.message);
        end
    end

    % plots
    figure()
    semilogy(ls_fac_vec, error(kk,:))
    xlabel('ls\_fac')
    title('Lengthscale-factor vs. rel. error')
    figure()
    plot(ls_fac_vec, loglike(kk,:))
    xlabel('ls\_fac')
    title('Lengthscale-factor vs. loglikelihood')
    
    % Best solution based on loglike
    [~, ii] = max(loglike(kk,:));
    ls_fac = ls_fac_vec(ii);
    disp(['Best ls_fac = ' num2str(ls_fac)])
    [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'kernel', kernel, 'maxit', maxit, 'abstol', abstol );
    disp(info)
    
    for ii=1:m
        figure()
        hold on
        plot(sol_ref.x, sol_ref.y(ii,:),'LineWidth', 1.5, 'LineStyle', '--', 'Color', color_vec(bvp5c_color,:))
        sd_fac = [8e4, 4e3];
        plot(sol.x, sol.y(ii,:)+sd_fac(kk)*sol.sd(ii,:),'LineWidth', 0.5,  'LineStyle', '-', 'Color', color_vec(end,:));
        plot(sol.x, sol.y(ii,:)-sd_fac(kk)*sol.sd(ii,:),'LineWidth', 0.5,  'LineStyle', '-', 'Color', color_vec(end,:));
        plot(sol.x, sol.y(ii,:),'LineWidth', 1.0,  'LineStyle', '-', 'Color', color_vec(our_color-1,:));

        xlabel('$t$', 'Interpreter','latex')
        ylabel(['$y_' num2str(ii) '(t)$'], 'Interpreter','latex')
        legend('hide')
        TikzFigure(['07_painleve_sol' num2str(kk) '_y' num2str(ii) '.tex'])
    end
    
    figure()
    for ii=1:m
        subplot(1,m,ii)
        plot(sol.x, sol.sd(ii,:), 'k');
        xlabel('t')
        ylabel(['sd(y_' num2str(ii-1) '(t))'])
    end
    
    figure()
    semilogy(sol.residual_rel)
    hold on
    semilogy(sol.step)
    
    
end





%%




