function [probs,funs,bcfuns,dfuns,dbcfuns,esolus,setoutputs,settolerancess]=bvpT34(ExtraArgs)
%   For a description of the input parameters see the report manual_testset
%
%   The structure of this code is based on the bvp testset, which can be
%   found here
%   https://archimede.dm.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=27
%  
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
         
%
probs    = @prob;
funs     = @fun;
bcfuns   = @bcfun;
dfuns    = @dfun;
dbcfuns  = @dbcfun;
esolus   = @esolu;
setoutputs = @setoutput;
settolerancess = @settolerances; 
if nargin==1 && ~isempty(ExtraArgs)
    if iscell(ExtraArgs)
       lambda=ExtraArgs{:};
    else
        lambda=ExtraArgs;
    end
end

    function [problm,type,m,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = prob()
        problm = 'bvpT34';
        type   = 'ODEBVP';
        t(1)   = 0;
        t(2)   = 10;
        m   = 2;     
        y0  = 0;        
        Linear  = 'off';
        numjac = 0;
        numbcjac = 0;
        Vectorized = 'on';
        JVectorized = 'on';
        solinit = bvpinit(linspace(t(1),t(2),11),repmat(y0,1,m));
    end

    function F = fun(X,Z,ExtraArgs)
        if nargin==3
            lambda=ExtraArgs;
        end
        F=zeros(size(Z));
        F(1,:) = Z(2,:);
        F(2,:) = (Z(1,:).*Z(1,:) - X)/lambda;
    end

    function  bc = bcfun(ya,yb,ExtraArgs)
        if nargin==3
            lambda=ExtraArgs;
        end
        C0 = [1,0;0,0];
        C1 = [0,0;1,0];
        Eta = [0;sqrt(10)];
        bc = C0*ya + C1*yb  - Eta;
    end

    function Df = dfun(X,Z,ExtraArgs)
        if nargin==3
            lambda=ExtraArgs;
        end
        ncomp=size(Z,1);
        nmsh=size(Z,2);
        Z=reshape(Z,ncomp,1,nmsh);
        Df=zeros(ncomp,ncomp,nmsh);
        Df(1,2,:) = 1.0e0;
        Df(2,1,:) = 2.e0*Z(1,1,:)/lambda;
    end

    function  [C0,C1] = dbcfun(ya,yb,ExtraArgs)
        if nargin==3
            lambda=ExtraArgs;
        end
        C0 = [1,0;0,0];
        C1 = [0,0;1,0];
    end

    function  tolvec = settolerances(tol)
     tolvec = tol;
    end 

    function [solref,printsolout,nindsol,indsol] = setoutput(plotsol)
      solref = 0; 
       if isempty(plotsol)
          printsolout = 0;
          nindsol = 1;
          indsol = 1;
       else
          printsolout = 1;
          nindsol = length(plotsol);
          indsol = plotsol;  
       end         
    end

    function Exact = esolu(X,ExtraArgs) 
        error('MATLAB:twpbvpc:bvp_examples:noexactsolution', 'No exact solution available');
    end
end