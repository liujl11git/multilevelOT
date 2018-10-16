function [m,phi] = W1PD_ML(dx, rho0, rho1, p, opts)
% An implement of the Algorithm 1M in the following paper
% J. Liu, W. Yin, W.C. Li, Y.T. Chow, "Multilevel Optimal Transport: 
% a Fast Approximation of Wasserstein-1 distances", submitted, 2018.
% 
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10
% 

%% paramters
if p ~= 1 && p ~= 2 && p~= 3, error('p should be 1,2, or 3.'); end

opts = checkopts(opts, size(rho0,1), dx);
mu=opts.mu;
tau = opts.tau;
level = opts.L;
verbose = opts.verbose;

%% initialization
rho0s = cell(1,level);
rho1s = cell(1,level);
rho0s{level} = rho0;
rho1s{level} = rho1;

normalizerho = @(x) x / sum(x(:));

for ll = (level-1):(-1):1
    h = dx * (2)^(level-ll);
    rho0s{ll} = restrict(rho0s{ll+1});
    rho1s{ll} = restrict(rho1s{ll+1});
    rho0s{ll} = normalizerho(rho0s{ll})/h/h;
    rho1s{ll} = normalizerho(rho1s{ll})/h/h;
end
Mx = size(rho0s{1},1);
My = size(rho0s{1},2);
m_next = zeros(Mx, My, 2);
phi_next = zeros(Mx, My);
ptimer = tic;
lasttime = 0;

%% main loop
for ll = 1:level
	h = dx * (2)^(level-ll);
	opts.mu = mu * (2)^(level-ll);
    opts.tau = tau * (2)^(level-ll);
    opts.m = m_next;
    opts.phi = phi_next;
    opts.verbose = 0;
	[m, phi, ~, iters] = W1PD(h, rho0s{ll}, rho1s{ll}, p, opts);
    timee = toc(ptimer) - lasttime;
    lasttime = toc(ptimer);
    
    if verbose == 1
        fprintf('l:%d\t iters:%d\t time:%f\n', ll, iters, timee);
    end
    
    if ll < level
        m_next = interpolate_m(m);
        phi_next = interpolate(phi);
    end
end
end


%% defalt options
function opts = checkopts(opts, M, h)

if isfield(opts, 'mu') == 0, opts.mu = h/2.7; end
if isfield(opts, 'tau') == 0, opts.tau = h/2.7; end
if isfield(opts, 'MaxIter') == 0, opts.MaxIter = 100000; end
if isfield(opts, 'tol') == 0, opts.tol = 1e-6; end
if isfield(opts, 'verbose') == 0, opts.verbose = 1; end
if isfield(opts, 'displayIter') == 0, opts.displayIter = 100; end
if isfield(opts, 'L') == 0
    L = log(M-1)/log(2) - 1;
    L = min(L,6);
    L = max(L,1);
    opts.L = L;
end

end

