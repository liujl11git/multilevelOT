function [m,gphi,potential] = W1PDHG_ML(dx, rho0, rho1, p,opts)
% An implement of the Algorithm 2M in the following paper
% J. Liu, W. Yin, W.C. Li, Y.T. Chow, "Multilevel Optimal Transport: 
% a Fast Approximation of Wasserstein-1 distances", submitted, 2018.
% 
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10
% 
% Note: This function depends on a mex file (emdPDHG.mexa64 for linux;
% emdPDHG.mexw64 for windows; emdPDHG.mexmaci64 for macOS)
% The mexfile is compiled with a C source code authored by 
% Matthew Jacobs @ UCLA (majaco@math.ucla.edu).

%% paramters
if p ~= 1 && p ~= 2 && p~= 3, error('p should be 1,2, or 3.'); end

opts = checkopts(opts, size(rho0,1));
tau = opts.tau;
mu = opts.mu;
MaxIter = opts.MaxIter;
tol = opts.tol;
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
u = zeros(Mx, My, 2);
gphi = zeros(Mx, My, 2);

%% main loop     
for ll = 1:level
       
    ptimer = tic;
	density = rho0s{ll}-rho1s{ll};
	[ux,uy,px,py,psix,psiy,potential,iters,ffttime, calctime] = ...
        emdPDHG(density, mu, tau, ...
        tol, MaxIter, u(:,:,1), u(:,:,2), gphi(:,:,1), gphi(:,:,2),p);
               
	u(:,:,1) = ux;
	u(:,:,2) = uy;
	gphi(:,:,1) = px;
	gphi(:,:,2) = py;
    
    % fprintf('fft-setup-time:%f  calc time:%f\n',ffttime, calctime);
    % the fft-setup time and calc-time are got by c++ clock()
    % which is less accurate than the MATLAB timer
    % because MATLAB timer uses platform specific functions
    
    if verbose == 1
        fprintf('l:%d\titers:%d\ttime:%f\n',ll, iters, toc(ptimer));
    end
    
    if ll < level
        u = interpolate_m2(u);
        gphi = interpolate_m2(gphi);
    end
end
    
u(:,:,1) = u(:,:,1) + psix;
u(:,:,2) = u(:,:,2) + psiy;
m = u;

end

%% default options
function opts = checkopts(opts, M)

if isfield(opts, 'mu') == 0, opts.mu = 1.0; end
if isfield(opts, 'tau') == 0, opts.tau = 1.0; end
if isfield(opts, 'MaxIter') == 0, opts.MaxIter = 100000; end
if isfield(opts, 'tol') == 0, opts.tol = 1e-5; end
if isfield(opts, 'verbose') == 0, opts.verbose = 1; end
if isfield(opts, 'L') == 0
    L = log(M-1)/log(2) - 1;
    L = min(L,6);
    L = max(L,1);
    opts.L = L;
end

end

%% interpolation
function [mf] = interpolate_m2(m)

Mx = size(m,1)-1; My = size(m,2)-1;
Mxf = Mx*2; Myf = My*2;
mf = zeros(Mxf+1, Myf+1, 2);

mf(1, 1:2:Myf+1, 1) = m(1, :, 1);
mf(2:2:Mxf, 1:2:Myf+1, 1) = m(2:Mx+1, :, 1);
mf(3:2:Mxf+1, 1:2:Myf+1, 1) = m(2:Mx+1, :, 1);
mf(1, 2:2:Myf, 1) = (m(1, 1:My, 1)+m(1, 2:My+1, 1))/2;
mf(2:2:Mxf, 2:2:Myf, 1) = (m(2:Mx+1, 1:My, 1)+m(2:Mx+1, 2:My+1, 1))/2;
mf(3:2:Mxf+1, 2:2:Myf, 1) = (m(2:Mx+1, 1:My, 1)+m(2:Mx+1, 2:My+1, 1))/2;

mf(1:2:Mxf+1, 1, 2) = m(:, 1, 2);
mf(1:2:Mxf+1, 2:2:Myf, 2) = m(:, 2:My+1, 2);
mf(1:2:Mxf+1, 3:2:Myf+1, 2) = m(:, 2:My+1, 2);
mf(2:2:Mxf, 1, 2) = (m(1:Mx, 1, 2)+m(2:Mx+1, 1, 2))/2;
mf(2:2:Mxf, 2:2:Myf, 2) = (m(1:Mx, 2:My+1, 2)+m(2:Mx+1, 2:My+1, 2))/2;
mf(2:2:Mxf, 3:2:Myf+1, 2) = (m(1:Mx, 2:My+1, 2)+m(2:Mx+1, 2:My+1, 2))/2;

end
