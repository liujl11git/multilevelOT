function [P] = ShrinkLinf(D, tau)
%%%%   min max(Px, Py) + 1/2/tau ((Px - Dx)^2 + (Py - Dy)^2)
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

P = D - ProjL1(D,tau);

end
