%  Purpose: Check influence of kernel lenghtscale for linear or nonlinear
%  problems 
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
ee = 1e-1;
[probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvpT19(ee);
[problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();

%%
M = 20;
npoints = 71;
kernel = @se_kernel; %@matern_52; % 
ls_fac_vec = logspace(log10(1.5),log10(15),M);
error = NaN(M,1);
loglike = NaN(M,1);

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-12); % default 'RelTol', 1e-3, 'AbsTol', 1e-6            
sol_ref = bvp4c(odefun,bcfun,solinit, options);
% solinit_ref = bvpinit(sol_ref.x,solinit.yinit);

for ii = 1:M
    ls_fac = ls_fac_vec(ii);
    try
        disp(['ls_fac = ' num2str(ls_fac)])
        [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'npoints', npoints, 'ls_fac', ls_fac, 'kernel', kernel);
        disp(info)
        try
            sol_ref_y = esolus(sol.x);
            temp = sol_ref_y(1,:);
        catch
            sol_ref_y = sol_ref.y;
            temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x);
        end
        error(ii) = norm(sol.y(1,:) - temp)/norm(temp);
        if 0 %abs(sol.loglike) > 1e4
            loglike(ii) = -Inf;
        else
            loglike(ii) = sol.loglike;
        end
        
    catch e %e is an MException struct
        fprintf(2,'The identifier was: %s\n',e.identifier);
        fprintf(2,'There was an error! The message was: %s\n',e.message);
    end
end
figure()
semilogy(ls_fac_vec, error)
xlabel('ls\_fac')
title('Lengthscale-factor vs. rel. error')
figure()
plot(ls_fac_vec, loglike)
xlabel('ls\_fac')
title('Lengthscale-factor vs. loglikelihood')

%% best solution
figure()
[~, ii] = min(error);
ls_fac = ls_fac_vec(ii);
disp(['Best ls_fac = ' num2str(ls_fac)])
[sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'npoints', npoints, 'ls_fac', ls_fac, 'kernel', kernel);
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
    plot(sol_ref_x, sol_ref_y(ii,:),'r*')
    hold on
    plot(sol.x, sol.y(ii,:), 'k-');
    hold on
    plot(sol.x, sol.y(ii,:)+sol.sd(ii,:), 'k:');
    hold on
    plot(sol.x, sol.y(ii,:)-sol.sd(ii,:), 'k:');
    xlabel('t')
    ylabel(['y_' num2str(ii-1) '(t)'])
end

figure()
for ii=1:m
    subplot(1,m,ii)
    plot(sol.x, sol.sd(ii,:), 'k');
    xlabel('t')
    ylabel(['sd(y_' num2str(ii-1) '(t))'])
end

%%

[~, ii] = max(loglike);
ls_fac = ls_fac_vec(ii);
disp(['Best ls_fac = ' num2str(ls_fac)])


