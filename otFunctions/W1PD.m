function [m, phi, optinf, iters] = W1PD(dx, rho0, rho1, p, opts)
% An implement of the algorithm 
% in the following paper (without parallelization)
% W. Li, E. Ryu, S. Osher, W. Yin, W. Gangbo, 
% " A Parallel Method for Earth's Movers Distance", 
% Journal of Scientific Computing, 2017.
% 
% Author: Jialin Liu (liujl11@math.ucla.edu) 
% Author: Wuchen Li (wcli@math.ucla.edu) 
% Modified: 2018-10-10

%% paramters
if p ~= 1 && p ~= 2 && p~= 3, error('p should be 1,2, or 3.'); end

mu = opts.mu; 
tau = opts.tau;
MaxIter = opts.MaxIter;
tol = opts.tol;
verbose = opts.verbose;
dispIter = opts.displayIter;

%% initialization
if(isfield(opts,'m') && isfield(opts,'phi'))
    m = opts.m;
    phi = opts.phi;
else
    Mx = size(rho0,1);
    My = size(rho0,2);
    m = zeros(Mx, My, 2); 
    phi=zeros(Mx, My);
end

info_now = [];
info_now.primal_fun = 0;
info_now.constraint = 0;
info_now.pd_residual = 0;
info_now.support=0;

outLoops = floor(MaxIter/dispIter);
optinf = repmat(info_now, outLoops, 1);

%% main loop
for k = 1:outLoops
    
     [m,phi,info_now,iters] = PDupdate(p, m, phi, dx, mu, tau, ...
                                rho0, rho1, dispIter,tol);
     
     if verbose == 1
            fprintf('k:%d\t constraint_res:%f pd_res:%e func:%f\n',...
            k*dispIter, info_now.constraint, ...
            info_now.pd_residual, info_now.primal_fun);
     end

     optinf(k) = info_now;
     
     if info_now.pd_residual < tol, break; end
    
end

iters = (k-1)*dispIter + iters;

return


%% single step update
function [newm,newphi,optinf,kk] = PDupdate(p,m,phi,dx,mu,tau, ...
                                    rho0,rho1,dispIters,tol)
newm = m;
newphi = phi;

for kk = 1:dispIters
    m = newm;
    phi = newphi;
    
    if p == 1, newm = ShrinkL1(m + mu*Gra(phi,dx), mu); end
    if p == 2, newm = ShrinkL2(m + mu*Gra(phi,dx), mu); end
    if p == 3, newm = ShrinkLinf(m + mu*Gra(phi,dx), mu); end
    divnewm = Div(newm,dx);
    divm = Div(m,dx);
    newphi = phi + tau*(rho1 - rho0 + 2*divnewm - divm);
    
    residual = 1/mu * ( (norm(m(:)-newm(:)))^2 );
    residual = residual + 1/tau * (norm(phi(:)-newphi(:))^2);
    residual = residual - 2 * sum(sum( (newphi - phi).*(divnewm - divm) ));
    residual = residual * dx * dx;

    if residual < tol, break; end
end

optinf = [];
if p == 1, optinf.primal_fun = PrimalFunL1(m, dx); end
if p == 2, optinf.primal_fun = PrimalFunL2(m, dx); end
if p == 3, optinf.primal_fun = PrimalFunLinf(m, dx); end
optinf.constraint = Constraint(m, rho0, rho1, dx);
optinf.pd_residual = residual;
optinf.support = nnz(m(:));

return
