function [mf] = interpolate_m(m)
% A function to calculate the interpolation of flux m.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

Mx = size(m,1)-1; My = size(m,2)-1;
Mxf = Mx*2; Myf = My*2;
mf = zeros(Mxf+1, Myf+1, 2);

mf(2:2:Mxf, 1:2:Myf+1, 1) = m(1:Mx, :, 1);
mf(1:2:Mxf-1, 1:2:Myf+1, 1) = m(1:Mx, :, 1);
mf(2:2:Mxf, 2:2:Myf, 1) = (m(1:Mx, 1:My, 1)+m(1:Mx, 2:My+1, 1))/2;
mf(1:2:Mxf-1, 2:2:Myf, 1) = (m(1:Mx, 1:My, 1)+m(1:Mx, 2:My+1, 1))/2;

mf(1:2:Mxf+1, 2:2:Myf, 2) = m(:, 1:My, 2);
mf(1:2:Mxf+1, 1:2:Myf-1, 2) = m(:, 1:My, 2);
mf(2:2:Mxf, 2:2:Myf, 2) = (m(1:Mx, 1:My, 2)+m(2:Mx+1, 1:My, 2))/2;
mf(2:2:Mxf, 1:2:Myf-1, 2) = (m(1:Mx, 1:My, 2)+m(2:Mx+1, 1:My, 2))/2;

end
