function [P] = ShrinkL2(D,tau)
%%%%   min sqrt(Px.^2 + Py.^2) + 1/2/tau ((Px - Dx)^2 + (Py - Dy)^2)
% Author: Jialin Liu (liujl11@math.ucla.edu) 
% Author: Wuchen Li (wcli@math.ucla.edu) 
% Modified: 2018-10-10

Dx = D(:,:,1);
Dy = D(:,:,2);
Dnorm = sqrt(Dx.^2 + Dy.^2);
idx = find(Dnorm > 0);
tmp = max(Dnorm - tau,0);

P = zeros([size(Dx,1), size(Dy,2), 2]);
Px = zeros(size(Dx));
Px(idx) = (Dx(idx)./Dnorm(idx)).*tmp(idx);
P(:,:,1) = Px;

Py = zeros(size(Dy));
Py(idx) = (Dy(idx)./Dnorm(idx)).*tmp(idx);
P(:,:,2) = Py;

end
