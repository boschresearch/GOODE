%  Purpose: Demonstrates the general process of solving nonlinear BVPs with
%  Jon's algorithm for linear problems

%  Author: Sc
%  Creation: 2018-07-19
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

close all;

%% problem 1: ey'' - y = 0

clear;
clf;
holdtime = 2.0;

DT = [0; 1];
y0 = 1;
y1 = 0;

e = 1e-2;

ys = @(t) (exp(-t/sqrt(e)) - exp((t-2)/sqrt(e)))/(1 - exp(-2/sqrt(e)));

n = 11;
N = 101;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

% ----- The kernel function
ker = se_kernel(1.0/(n-1)*1.0, 1.);

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT), kxy(ker, DT, T), kxddy(ker, DT, T);
     kxy(ker, T, DT), kxy(ker, T, T), kxddy(ker, T, T);
     kddxy(ker, T, DT), kddxy(ker, T, T), kddxddy(ker, T, T)];
H = [1, 0, zeros(1, n);
     0, 1, zeros(1, n);
     zeros(n,2), diag(-1*ones(n,1));
     zeros(n,2), diag(e*ones(n,1))]';
   
Y = [y0; y1; zeros(n,1)];

k = [kxy(ker, t, DT), kxy(ker, t, T), kxddy(ker, t, T)];

L = chol(H * G * H' + 1e-10*eye(size(H*H')));

mu = k * H' * (L \ (L' \ Y));
plot(t, ys(t));  hold on; plot(t, mu); pause(holdtime);

k2 = [kxy(ker, DT, t); kxy(ker, T, t); kddxy(ker, T, t)];
k_post = kxy(ker, t, t) - k * H' * (L \ (L' \ (H* k2)));
plot(t, mu+sqrt(diag(k_post)), t,mu-sqrt(diag(k_post)) )

cond(G)
cond(H * G * H')
cond(H)

%% problem 2: ey'' - y' = 0

clear;
clf;
holdtime = 2.0;

DT = [0; 1];
y0 = 1;
y1 = 0;

e = 7e-3;

ys = @(t) (1 - exp((t-1)/e))/(1-exp(-1/e));

n = 71;
N = 51;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

% ----- The kernel function
ker = se_kernel(0.05, 1.);

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT), kxdy(ker, DT, T), kxddy(ker, DT, T);
     kdxy(ker, T, DT), kdxdy(ker, T, T), kdxddy(ker, T, T);
     kddxy(ker, T, DT), kddxdy(ker, T, T), kddxddy(ker, T, T)];
H = [1, 0, zeros(1, n);
     0, 1, zeros(1, n);
     zeros(n,2), diag(-1*ones(n,1));
     zeros(n,2), diag(e*ones(n,1))]';
   
Y = [y0; y1; zeros(n,1)];

k = [kxy(ker, t, DT), kxdy(ker, t, T), kxddy(ker, t, T)];

L = chol(H * G * H' + 1e-10*eye(size(H*H')));

mu = k * H' * (L \ (L' \ Y));
plot(t, ys(t));  hold on; plot(t, mu); pause(holdtime);

%% problem 4: ey'' + y' - (1+e)y= 0

clear;
clf;
holdtime = 2.0;

e = 5e-2;

DT = [-1; 1];
y0 = 1 + exp(-2);
y1 = 1+exp(-2*(1+e)/e);

ys = @(t) exp(t-1)+exp(-(1+e)*(1+t)/e);

n = 51;
N = 51;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

% ----- The kernel function
ker = se_kernel(0.05, 1.);

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT), kxy(ker, DT, T), kxdy(ker, DT, T), kxddy(ker, DT, T);
     kxy(ker, T, DT), kxy(ker, T, T), kxdy(ker, T, T), kxddy(ker, T, T);
     kdxy(ker, T, DT), kdxy(ker, T, T), kdxdy(ker, T, T), kdxddy(ker, T, T);
     kddxy(ker, T, DT), kddxy(ker, T, T), kddxdy(ker, T, T), kddxddy(ker, T, T)];
H = [1, 0, zeros(1, n);
     0, 1, zeros(1, n);
     zeros(n,2), diag(-(1+e)*ones(n,1));
     zeros(n,2), diag(1*ones(n,1));
     zeros(n,2), diag(e*ones(n,1))]';
   
Y = [y0; y1; zeros(n,1)];

k = [kxy(ker, t, DT), kxy(ker, t, T), kxdy(ker, t, T), kxddy(ker, t, T)];

L = chol(H * G * H' + 1e-10*eye(size(H*H')));

mu = k * H' * (L \ (L' \ Y));
plot(t, ys(t));  hold on; plot(t, mu); pause(holdtime);
  
%% problem 6: ey'' + ty' = -e pi^2cos(pi t) - pi t sin(pi t)

clear;
clf;
holdtime = 2.0;

DT = [-1; 1];
y0 = -2;
y1 = 0;

e = 1e-3;

ys = @(t) cos(pi*t) + erf(t/sqrt(2*e))/erf(1/sqrt(2*e));

n = 71;
N = 101;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

% ----- The kernel function
ker = se_kernel(0.05, 1.);

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT), kxdy(ker, DT, T), kxddy(ker, DT, T);
     kdxy(ker, T, DT), kdxdy(ker, T, T), kdxddy(ker, T, T);
     kddxy(ker, T, DT), kddxdy(ker, T, T), kddxddy(ker, T, T)];
H = [1, 0, zeros(1, n);
     0, 1, zeros(1, n);
     zeros(n,2), diag(T);
     zeros(n,2), diag(e*ones(n,1))]';
   
Y = [y0; y1; -e*pi^2*cos(pi*T) - pi*T.*sin(pi*T)];

k = [kxy(ker, t, DT), kxdy(ker, t, T), kxddy(ker, t, T)];

L = chol(H * G * H' + 1e-10*eye(size(H*H')));

mu = k * H' * (L \ (L' \ Y));
plot(t, ys(t));  hold on; plot(t, mu); pause(holdtime);

%% Nonlinear problems

%                               $$\ $$\                                         
%                               $$ |\__|                                        
% $$$$$$$\   $$$$$$\  $$$$$$$\  $$ |$$\ $$$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\  
% $$  __$$\ $$  __$$\ $$  __$$\ $$ |$$ |$$  __$$\ $$  __$$\  \____$$\ $$  __$$\ 
% $$ |  $$ |$$ /  $$ |$$ |  $$ |$$ |$$ |$$ |  $$ |$$$$$$$$ | $$$$$$$ |$$ |  \__|
% $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |$$ |$$ |  $$ |$$   ____|$$  __$$ |$$ |      
% $$ |  $$ |\$$$$$$  |$$ |  $$ |$$ |$$ |$$ |  $$ |\$$$$$$$\ \$$$$$$$ |$$ |      
% \__|  \__| \______/ \__|  \__|\__|\__|\__|  \__| \_______| \_______|\__|                                                                                    

%% problem 19: ey'' + exp(y)y' - .5*pi sin(pi t/2)exp(2y) = 0

clear;
clf;
itertime = 1.0;
holdtime = 2.0;

e = 1e-2; % minimum value: 9e-3

DT = [0; 1];
y0 = 0;
y1 = 0;

% ys = @(t) 1 + e*log(cosh((t-0.745)/e));

n = 51;
N = 101;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

f   = @(t,y,dy) (0.5*pi*sin(0.5*pi*t(:)).*exp(2*y(:)) - exp(y(:)).*dy(:))/e;
fy  = @(t,y,dy) exp(y).*(pi*sin(0.5*pi*t(:)).*exp(y(:)) - dy(:))/e;
fdy = @(t,y,dy) -exp(y(:))/e;

% ----- The mean functions
m = @(t) y0 + (y1 - y0) * (t - DT(1)) ./ (DT(2) - DT(1));
dm = @(t) ones(size(t, 1), 1) * (y1 - y0)' ./ (DT(2) - DT(1)); % T x D

% ----- the true solution

odeinit = struct();
odeinit.x = T';
odeinit.y = [m(T)'; dm(T)'];

odefun = @(t,y) [y(2,:); f(t,y(1,:),y(2,:))];
bcfun  = @(ya,yb) [ya(1); yb(1)];

sol = tom(odefun, bcfun, odeinit);

% ----- The kernel function
ker = se_kernel(0.05, 1.);

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT),  kxy(ker, DT, T),  kxdy(ker, DT, T),  kxddy(ker, DT, T);
     kxy(ker, T, DT),   kxy(ker, T, T),   kxdy(ker, T, T),   kxddy(ker, T, T);
     kdxy(ker, T, DT),  kdxy(ker, T, T),  kdxdy(ker, T, T),  kdxddy(ker, T, T);
     kddxy(ker, T, DT), kddxy(ker, T, T), kddxdy(ker, T, T), kddxddy(ker, T, T)];

kT = [kxy(ker, T, DT),  kxy(ker, T, T),  kxdy(ker, T, T),  kxddy(ker, T, T);
      kdxy(ker, T, DT), kdxy(ker, T, T), kdxdy(ker, T, T), kdxddy(ker, T, T)];

k  = [kxy(ker, t, DT), kxy(ker, t, T), kxdy(ker, t, T), kxddy(ker, t, T)];

% ----- iteration variables

Y0 = [y0; y1; f(T,m(T),dm(T)) - fy(T,m(T),dm(T)).*m(T) - fdy(T,m(T),dm(T)).*dm(T)];

H0 = [1, 0, zeros(1, n);
      0, 1, zeros(1, n);
      zeros(n,2), diag(-fy(T,m(T),dm(T)));
      zeros(n,2), diag(-fdy(T,m(T),dm(T)));
      zeros(n,2), diag(ones(n,1))]';

L0 = chol(H0 * G * H0' + 1e-10*eye(size(H0,1)));

ai = H0' * (L0 \ (L0' \ Y0));

muTi = reshape(kT * ai,n,2);

plot(sol.x, sol.y(1,:)); hold on; plot(t, k * ai); pause(itertime);

for i=1:5
  Yi = [y0; y1; f(T,muTi(:,1),muTi(:,2)) ...
                - fy(T,muTi(:,1),muTi(:,2)).*muTi(:,1) ...
                - fdy(T,muTi(:,1),muTi(:,2)).*muTi(:,2)];
  Hi = [1, 0, zeros(1, n);
        0, 1, zeros(1, n);
        zeros(n,2), diag(- fy(T,muTi(:,1),muTi(:,2)));
        zeros(n,2), diag(- fdy(T,muTi(:,1),muTi(:,2)));
        zeros(n,2), diag(ones(n,1))]';
      
  Li = chol(Hi * G * Hi' + 1e-10*eye(size(Hi,1)));
  ai = Hi' * (Li \ (Li' \ Yi));
  
  muTi = reshape(kT * ai,n,2);
  
  plot(t, k * ai); pause(itertime);
end

pause(holdtime-itertime);

%% problem 20: ey'' + (y')^2 = 1

clear;
clf;
itertime = .5;
holdtime = 2.0;

e = 5e-2; % minimum value 4e-2

DT = [0; 1];
y0 = 1 + e*log(cosh(-0.745/e));
y1 = 1 + e*log(cosh(0.255/e));

ys = @(t) 1 + e*log(cosh((t-0.745)/e));

n = 61;
N = 101;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

% ----- The mean functions
m = @(t) y0 + (y1 - y0) * (t - DT(1)) ./ (DT(2) - DT(1));
dm = @(t) ones(size(t, 1), 1) * (y1 - y0)' ./ (DT(2) - DT(1)); % T x D

% ----- The kernel function
ker = se_kernel(0.05, 1.);

f   = @(t,dy) (1-dy(:).^2)/e;
fdy = @(t,dy) -2*dy(:)/e;

%  ----- The kernel parts for the GP
G  = [kxy(ker, DT, DT),  kxdy(ker, DT, T),  kxddy(ker, DT, T);
      kdxy(ker, T, DT),  kdxdy(ker, T, T),  kdxddy(ker, T, T);
      kddxy(ker, T, DT), kddxdy(ker, T, T), kddxddy(ker, T, T)];
kT = [kdxy(ker, T, DT),  kdxdy(ker, T, T),  kdxddy(ker, T, T)];

k  = [kxy(ker, t, DT), kxdy(ker, t, T), kxddy(ker, t, T)];

% ----- iteration variables

Y0 = [y0; y1; f(T,dm(T)) - fdy(T,dm(T)).*dm(T)];

H0 = [1, 0, zeros(1, n);
      0, 1, zeros(1, n);
      zeros(n,2), diag(-fdy(T,dm(T)));
      zeros(n,2), diag(ones(n,1))]';
L0 = chol(H0 * G * H0' + 1e-10*eye(size(H0,1)));

ai = H0' * (L0 \ (L0' \ Y0));

muTi = kT * ai;

plot(t, ys(t)); hold on; plot(t, k * ai); pause(itertime);

for i=1:10
  Yi = [y0; y1;f(T,muTi) - fdy(T,muTi).*muTi];
  Hi = [1, 0, zeros(1, n);
        0, 1, zeros(1, n);
        zeros(n,2), diag(-fdy(T,muTi));
        zeros(n,2), diag(ones(n,1))]';
      
  Li = chol(Hi * G * Hi' + 1e-10*eye(size(Hi,1)));
  ai = Hi' * (Li \ (Li' \ Yi));
  
  muTi = kT * ai;
  
  plot(t, k * ai); pause(itertime);
end

pause(holdtime-itertime);

%% problem 21: ey'' = y + y^2 - exp(-2 t e^-.5)

clear;
clf;
itertime = 1.;
holdtime = 2.0;

e = 1e-4; % minimum value 5e-5

DT = [0; 1];
y0 = 1;
y1 = exp(-1/sqrt(e));

ys = @(t) exp(-t/sqrt(e));

n = 51;
N = 101;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

% ----- The mean functions
m = @(t) y0 + (y1 - y0) * (t - DT(1)) ./ (DT(2) - DT(1));

% ----- The kernel function
ker = se_kernel(0.05, 1.);

f  = @(t,y) y(:)/e + y(:).^2/e - exp(-2*t(:)/sqrt(e))/e;
fy = @(t,y) 1/e + 2*y(:)/e;

%  ----- The kernel parts for the GP
G  = [kxy(ker, DT, DT),  kxy(ker, DT, T),  kxddy(ker, DT, T);
      kxy(ker, T, DT),   kxy(ker, T, T),   kxddy(ker, T, T);
      kddxy(ker, T, DT), kddxy(ker, T, T), kddxddy(ker, T, T)];
kT = [kxy(ker, T, DT), kxy(ker, T, T), kxddy(ker, T, T)];

k  = [kxy(ker, t, DT), kxy(ker, t, T), kxddy(ker, t, T)];

% ----- iteration variables

Y0 = [y0; y1; f(T,m(T)) - fy(T,m(T)).*m(T)];

H0 = [1, 0, zeros(1, n);
      0, 1, zeros(1, n);
      zeros(n,2), diag(-fy(T,m(T)));
      zeros(n,2), diag(ones(n,1))]';
L0 = chol(H0 * G * H0' + 1e-10*eye(size(H0,1)));

ai = H0' * (L0 \ (L0' \ Y0));

muTi = kT * ai;

plot(t, ys(t)); hold on; plot(t, k * ai); pause(itertime);

for i=1:5
  Yi = [y0; y1; f(T,muTi) - fy(T,muTi).*muTi];
  Hi = [1, 0, zeros(1, n);
        0, 1, zeros(1, n);
        zeros(n,2), diag(-fy(T,muTi));
        zeros(n,2), diag(ones(n,1))]';
      
  Li = chol(Hi * G * Hi' + 1e-10*eye(size(Hi,1)));
  ai = Hi' * (Li \ (Li' \ Yi));
  
  muTi = kT * ai;
  
  plot(t, k * ai); pause(itertime);
end

pause(holdtime-itertime);
 
%% problem 22: ey'' + y' + y^2 = 0
 
clear;
clf;
itertime = 1.0;
holdtime = 2.0;

e = 7e-3; % minimum value ~7e-3

DT = [0; 1];
y0 = 0;
y1 = 0.5;

% ys = @(t) 1 + e*log(cosh((t-0.745)/e));

n = 61;
N = 101;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

f   = @(t,y,dy) (-(y(:).^2) - dy(:))/e;
fy  = @(t,y,dy) -2*y(:)/e;
fdy = @(t,y,dy) -ones(size(dy(:)))/e;

% ----- The mean functions
m = @(t) y0 + (y1 - y0) * (t - DT(1)) ./ (DT(2) - DT(1));
dm = @(t) ones(size(t, 1), 1) * (y1 - y0)' ./ (DT(2) - DT(1)); % T x D

% ----- the true solution

odeinit = struct();
odeinit.x = T';
odeinit.y = [m(T)'; dm(T)'];

odefun = @(t,y) [y(2,:); f(t,y(1,:),y(2,:))];
bcfun  = @(ya,yb) [ya(1) - y0; yb(1) - y1];

sol = tom(odefun, bcfun, odeinit);

% ----- The kernel function
ker = se_kernel(0.05, 1.);

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT),  kxy(ker, DT, T),  kxdy(ker, DT, T),  kxddy(ker, DT, T);
     kxy(ker, T, DT),   kxy(ker, T, T),   kxdy(ker, T, T),   kxddy(ker, T, T);
     kdxy(ker, T, DT),  kdxy(ker, T, T),  kdxdy(ker, T, T),  kdxddy(ker, T, T);
     kddxy(ker, T, DT), kddxy(ker, T, T), kddxdy(ker, T, T), kddxddy(ker, T, T)];

kT = [kxy(ker, T, DT),  kxy(ker, T, T),  kxdy(ker, T, T),  kxddy(ker, T, T);
      kdxy(ker, T, DT), kdxy(ker, T, T), kdxdy(ker, T, T), kdxddy(ker, T, T)];

k  = [kxy(ker, t, DT), kxy(ker, t, T), kxdy(ker, t, T), kxddy(ker, t, T)];

% ----- iteration variables

Y0 = [y0; y1; f(T,m(T),dm(T)) - fy(T,m(T),dm(T)).*m(T) - fdy(T,m(T),dm(T)).*dm(T)];

H0 = [1, 0, zeros(1, n);
      0, 1, zeros(1, n);
      zeros(n,2), diag(-fy(T,m(T),dm(T)));
      zeros(n,2), diag(-fdy(T,m(T),dm(T)));
      zeros(n,2), diag(ones(n,1))]';

L0 = chol(H0 * G * H0' + 1e-10*eye(size(H0,1)));

ai = H0' * (L0 \ (L0' \ Y0));

muTi = reshape(kT * ai,n,2);

plot(sol.x, sol.y(1,:)); hold on; plot(t, k * ai); pause(itertime);

for i=1:5
  Yi = [y0; y1; f(T,muTi(:,1),muTi(:,2)) ...
                - fy(T,muTi(:,1),muTi(:,2)).*muTi(:,1) ...
                - fdy(T,muTi(:,1),muTi(:,2)).*muTi(:,2)];
  Hi = [1, 0, zeros(1, n);
        0, 1, zeros(1, n);
        zeros(n,2), diag(- fy(T,muTi(:,1),muTi(:,2)));
        zeros(n,2), diag(- fdy(T,muTi(:,1),muTi(:,2)));
        zeros(n,2), diag(ones(n,1))]';
      
  Li = chol(Hi * G * Hi' + 1e-10*eye(size(Hi,1)));
  ai = Hi' * (Li \ (Li' \ Yi));
  
  muTi = reshape(kT * ai,n,2);
  
  plot(t, k * ai); pause(itertime);
end

pause(holdtime-itertime);

%% problem 25: ey'' + yy' - y = 0

clear;
clf;
itertime = 1.0;
holdtime = 2.0;

e = 5e-4; % minimum value ~1e-5

DT = [0; 1];
y0 = -1/3;
y1 = 1/3;

% ys = @(t) 1 + e*log(cosh((t-0.745)/e));

n = 41;
N = 101;
T = linspace(DT(1), DT(2), n)';
t = linspace(DT(1), DT(2), N)';

f   = @(t,y,dy) (y(:)-y(:).*dy(:))/e;
fy  = @(t,y,dy) (1 - dy(:))/e;
fdy = @(t,y,dy) -y(:)/e;

% ----- The mean functions
m = @(t) y0 + (y1 - y0) * (t - DT(1)) ./ (DT(2) - DT(1));
dm = @(t) ones(size(t, 1), 1) * (y1 - y0)' ./ (DT(2) - DT(1)); % T x D

% ----- the true solution

odeinit = struct();
odeinit.x = T';
odeinit.y = [m(T)'; dm(T)'];

odefun = @(t,y) [y(2,:); f(t,y(1,:),y(2,:))];
bcfun  = @(ya,yb) [ya(1) - y0; yb(1) - y1];

sol = tom(odefun, bcfun, odeinit);

% ----- The kernel function
ker = se_kernel(0.05, 1.);

%  ----- The kernel parts for the GP
G = [kxy(ker, DT, DT),  kxy(ker, DT, T),  kxdy(ker, DT, T),  kxddy(ker, DT, T);
     kxy(ker, T, DT),   kxy(ker, T, T),   kxdy(ker, T, T),   kxddy(ker, T, T);
     kdxy(ker, T, DT),  kdxy(ker, T, T),  kdxdy(ker, T, T),  kdxddy(ker, T, T);
     kddxy(ker, T, DT), kddxy(ker, T, T), kddxdy(ker, T, T), kddxddy(ker, T, T)];

kT = [kxy(ker, T, DT),  kxy(ker, T, T),  kxdy(ker, T, T),  kxddy(ker, T, T);
      kdxy(ker, T, DT), kdxy(ker, T, T), kdxdy(ker, T, T), kdxddy(ker, T, T)];

k  = [kxy(ker, t, DT), kxy(ker, t, T), kxdy(ker, t, T), kxddy(ker, t, T)];

% ----- iteration variables

Y0 = [y0; y1; f(T,m(T),dm(T)) - fy(T,m(T),dm(T)).*m(T) - fdy(T,m(T),dm(T)).*dm(T)];

H0 = [1, 0, zeros(1, n);
      0, 1, zeros(1, n);
      zeros(n,2), diag(-fy(T,m(T),dm(T)));
      zeros(n,2), diag(-fdy(T,m(T),dm(T)));
      zeros(n,2), diag(ones(n,1))]';

L0 = chol(H0 * G * H0' + 1e-10*eye(size(H0,1)));

ai = H0' * (L0 \ (L0' \ Y0));

muTi = reshape(kT * ai,n,2);

plot(sol.x, sol.y(1,:)); hold on; plot(t, k * ai); pause(itertime);

for i=1:5
  Yi = [y0; y1; f(T,muTi(:,1),muTi(:,2)) ...
                - fy(T,muTi(:,1),muTi(:,2)).*muTi(:,1) ...
                - fdy(T,muTi(:,1),muTi(:,2)).*muTi(:,2)];
  Hi = [1, 0, zeros(1, n);
        0, 1, zeros(1, n);
        zeros(n,2), diag(- fy(T,muTi(:,1),muTi(:,2)));
        zeros(n,2), diag(- fdy(T,muTi(:,1),muTi(:,2)));
        zeros(n,2), diag(ones(n,1))]';
      
  Li = chol(Hi * G * Hi' + 1e-10*eye(size(Hi,1)));
  ai = Hi' * (Li \ (Li' \ Yi));
  
  muTi = reshape(kT * ai,n,2);
  
  plot(t, k * ai); pause(itertime);
end

pause(holdtime-itertime);
