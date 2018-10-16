function [P] = ProjL1(D, tau)
%%%%   projection to \|D\|_1 <= tau
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

proj_pos_x = @(Dx,Dy,tau) (Dx-Dy+tau)/2;
proj_pos_y = @(Dx,Dy,tau) (Dy-Dx+tau)/2;

Dx = D(:,:,1);
Dy = D(:,:,2);

% phases
ind0 = Dy+Dx <= tau & Dy+Dx > -tau & Dy-Dx <= tau & Dy-Dx > -tau;
ind1 = Dy+Dx >  tau &                Dy-Dx <= tau & Dy-Dx > -tau;
ind2 = Dy+Dx <= tau & Dy+Dx > -tau & Dy-Dx >  tau;
ind3 = Dy+Dx < -tau &                Dy-Dx <= tau & Dy-Dx > -tau;
ind4 = Dy+Dx <= tau & Dy+Dx > -tau & Dy-Dx < -tau;
ind5 = Dy+Dx >  tau &                Dy-Dx >  tau;
ind6 = Dy+Dx < -tau &                Dy-Dx >  tau;
ind7 = Dy+Dx < -tau &                Dy-Dx < -tau;
ind8 = Dy+Dx >  tau &                Dy-Dx < -tau;

Px = zeros(size(Dx));
Py = Px;

Px(ind1) = proj_pos_x(Dx(ind1),Dy(ind1),tau);
Py(ind1) = proj_pos_y(Dx(ind1),Dy(ind1),tau);

Px(ind2) = -proj_pos_x(-Dx(ind2),Dy(ind2),tau);
Py(ind2) = proj_pos_y(-Dx(ind2),Dy(ind2),tau);

Px(ind3) = -proj_pos_x(-Dx(ind3),-Dy(ind3),tau);
Py(ind3) = -proj_pos_y(-Dx(ind3),-Dy(ind3),tau);

Px(ind4) = proj_pos_x(Dx(ind4),-Dy(ind4),tau);
Py(ind4) = -proj_pos_y(Dx(ind4),-Dy(ind4),tau);

Px(ind0) = Dx(ind0);
Py(ind0) = Dy(ind0);

Px(ind5) = 0;
Py(ind5) = tau;

Px(ind6) = -tau;
Py(ind6) = 0;

Px(ind7) = 0;
Py(ind7) = -tau;

Px(ind8) = tau;
Py(ind8) = 0;

P = zeros(size(D));
P(:,:,1) = Px;
P(:,:,2) = Py;

end
