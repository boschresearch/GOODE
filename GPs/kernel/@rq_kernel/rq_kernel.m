function k = rq_kernel (alpha, lambda)
% RQ_KERNEL - Represents a rational quadratic kernel
%
% Inputs:
%   alpha  - length-scale mixture
%   lambda - length-scale
%
% Returns:
%   k - Kernel object

% rq_kernel.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-04-28
% Version: 0.1

  k = struct ();
  
  if nargin < 2
    lambda = 1.;
    
    if nargin < 1
      alpha = 1.;
    end
  end
  
  k.alpha  = alpha;
  k.lambda = lambda;
  
  k = class(k, 'rq_kernel');

end