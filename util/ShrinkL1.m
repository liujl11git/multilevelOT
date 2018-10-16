function [P] = ShrinkL1(D, tau)
%%%%   min |Px| + |Py| + 1/2/tau ((Px - Dx)^2 + (Py - Dy)^2)
% Author: Jialin Liu (liujl11@math.ucla.edu) 
% Author: Wuchen Li (wcli@math.ucla.edu) 
% Modified: 2018-10-10

Dx = D(:,:,1);
Dy = D(:,:,2);
Px = sign(Dx).*max(abs(Dx) - tau,0);
Py = sign(Dy).*max(abs(Dy) - tau,0);
P(:,:,1) = Px;
P(:,:,2) = Py;

end
