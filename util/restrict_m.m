function mc = restrict_m(m)
% A function to calculate the downsample of flux m.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

Mx = size(m,1)-1; My = size(m,2)-1;
Mxc = Mx/2; Myc = My/2;
mc = zeros(Mxc+1, Myc+1, 2);

mc(1:Mxc, :, 1) = (m(1:2:(Mx-1),1:2:My+1,1) + m(2:2:(Mx),1:2:My+1,1))/2;
mc(:, 1:Myc, 2) = (m(1:2:Mx+1, 1:2:Mx-1, 2) + m(1:2:Mx+1, 2:2:Mx, 2))/2;

end
