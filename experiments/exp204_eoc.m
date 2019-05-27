%  Purpose: EOC Plots of GP BVP solver
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

problems_index = [1, 2, 20, 21]; % currently just for one single problem (at least plots)

nbvp = numel(problems_index);
for ii = 1:nbvp
    bvp_list{ii} = str2func(['bvpT' int2str(problems_index(ii))]);
end

% hyperparameters
kernel = @se_kernel; %@matern_52; % 
npoints_vec = [5,9,13,19,27,41,51,71,91,111]; 
N = length(npoints_vec);
ee = 1e-1;
M = 40;
%ls_fac_vec = linspace(1.5,15,M);
ls_fac_vec = logspace(0,log10(15),M);



error = NaN(nbvp,M,N);
loglike = NaN(nbvp,M,N);
sigma2 = NaN(nbvp,M,N);

options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-6);
solver = @bvp5c;

esolus_mask = zeros(nbvp,1);


for ii = 1:nbvp
    disp(['problem ' num2str(problems_index(ii))])
    [probs,odefun,bcfun,dodefun,dbcfun,esolus,setoutputs,settolerancess]=bvp_list{ii}(ee);
    [problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = probs();
    
    % compute reference
    sol_ref = solver(odefun,bcfun,solinit, options);
    
    % iterate over npoints_vec
    for kk = 1:N
        npoints = npoints_vec(kk);
        
        % iterate over ls_fac_vec
        for jj = 1:M
            ls_fac = ls_fac_vec(jj);
            try
                %disp(['ls_fac = ' num2str(ls_fac)])
                [sol, info] = gpbvp_nonlinear_2(odefun, bcfun, dodefun, dbcfun, solinit, 'ls_fac', ls_fac, 'npoints', npoints, 'kernel', kernel);
                %disp(info)
                try % exact solution
                    sol_ref_y = esolus(sol.x);
                    temp = sol_ref_y(1,:);
                    esolus_mask(ii) = 1;
                catch e % computed reference solution
                    sol_ref_y = sol_ref.y;
                    temp = interp1(sol_ref.x,sol_ref.y(1,:),sol.x);
                    fprintf(2,'There was an error! The message was: %s\n',e.message);
                end
                error(ii,jj,kk) = norm(sol.y(1,:) - temp)/norm(temp);
                sigma2(ii,jj,kk) = mean(sol.sd(1,:));
                % check if loglike tends to -Inf
                loglike_max = max(loglike(ii,1:jj-1));
                rel_change = abs(loglike_max-sol.loglike)/abs(loglike_max);
                if rel_change> 1e4
                    disp(rel_change)
                    break
                else
                    loglike(ii,jj,kk) = sol.loglike;
                end
            catch e %e is an MException struct
                fprintf(2,'The identifier was: %s\n',e.identifier);
                fprintf(2,'There was an error! The message was: %s\n',e.message);
            end
        end
    end
end


%% plots
close all;

linestyle = {':', ':', ':', ':', ':'};
our_color = 7;

for ii = 1:nbvp
    
    % Select best according to error
    [error_min, ind] = min(error(ii,:,:));
    linear_ind = sub2ind(size(error),ii*ones(1,N), reshape(ind,1,N), 1:N);
    sigma2_min = sigma2(linear_ind);
    error_min = reshape(error_min,1,N);
    
    figure()
    loglog(npoints_vec, error_min, 'o', 'LineStyle', linestyle{ii}, 'LineWidth', 1.5, 'Color', color_vec(our_color,:), 'DisplayName', ['problem ' num2str(problems_index(ii))])
    xlabel('N')
    ylabel('RelErr')
    
    % Select best according to loglike
    [~, ind] = max(loglike(ii,:,:));
    linear_ind = sub2ind(size(error),ii*ones(1,N), reshape(ind,1,N), 1:N);
    error_opt_loglike = error(linear_ind);
    hold on
    loglog(npoints_vec, error_opt_loglike, '^', 'LineStyle', linestyle{ii}, 'LineWidth', 1.5, 'Color', color_vec(our_color-1,:), 'DisplayName', ['problem ' num2str(problems_index(ii))])
    % legend('error', 'error loglike')
    
    TikzFigure(['05_convergence_prob_' num2str(problems_index(ii)) '.tex'])
   

    figure()
    semilogy(npoints_vec, sigma2_min, 'o', 'LineStyle', linestyle{ii}, 'Color', color_vec(our_color,:), 'DisplayName', ['problem ' num2str(problems_index(ii))])
    xlabel('N')
    ylabel('Mean(sigma)')
    % title(['problem ' num2str(problems_index(ii))])
    
    % calculate eoc % err = C*h^p
    % log(err) = log(C) + p * log(h)
    % log(err1/err2) = p * log(h1/h2)
%     T = solinit.x;
%     h = (T(end) - T(1))./ (npoints_vec-1);
%     eoc = log(error_min(1:end-1)./error_min(2:end)) ./ log(h(1:end-1)./h(2:end));
%     figure(1)
%     subplot(2,2,4)
%     plot(npoints_vec(2:end),eoc)
%     title(['eoc problem ' num2str(problems_index(ii))])
    
    % plot lengthscale
%     figure(2)
%     subplot(2,2,3)
%     plot(npoints_vec, h.*ls_fac_vec(reshape(ind,1,N)))
%     xlabel('N')
%     ylabel('length scale')
%     title(['problem ' num2str(problems_index(ii))])




end