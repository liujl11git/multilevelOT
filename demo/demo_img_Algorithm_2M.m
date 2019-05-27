% A demo of W1PDHG_ML (Wasserstein-1 Prox-PDHG MultiLevel), 
% which is an implementation of Algorithm 2M in the following paper
% J. Liu, W. Yin, W.C. Li, Y.T. Chow, "Multilevel Optimal Transport: 
% a Fast Approximation of Wasserstein-1 distances", submitted, 2018.
% 
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10
% 

%% problem settings
[h,rho0,rho1,x,y] = ReadIMG('Cameraman.png','Lake.png');

p = 2; %p means the "groud metric" in the paper. 
% p = 1,2 or 3. p=3 means p=inf

%% algorthim parameters
opts = [];
opts.tol = 1e-5; % tolerance for fixed-point-residual
opts.verbose = 1; % display metrics
opts.L = 6; % number of Levels

%% calculation
tic;
[m,graphi,phi] = W1PDHG_ML(h, rho0, rho1, p, opts);
toc;

%% displaying
fprintf('energy:%f\n',PrimalFunL2(m, h));

PlotFlow(x,y,m,rho0-rho1);
PlotPotential(x,y,phi);
