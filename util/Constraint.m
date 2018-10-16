function [c]=Constraint(m,rho0, rho1, dx)
% A function to calculate the constraint residual
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

res = rho1-rho0+Div(m,dx);
c = norm(res(:))/norm(rho1(:)-rho0(:));
end
