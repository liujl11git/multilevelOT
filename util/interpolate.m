function vf = interpolate(v)
% A function to calculate the interpolation of potential v.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

Mx = size(v,1)-1;  My = size(v,2)-1;
Mxf = Mx*2; Myf = My*2;
vf = zeros(Mxf+1, Myf+1);

vf(1:2:Mxf+1,1:2:Myf+1) = v;
vf(1:2:Mxf+1,2:2:Myf) = (v(:,1:My)+v(:,2:My+1))/2;
vf(2:2:Mxf, 1:2:Myf+1) = (v(1:Mx,:)+v(2:Mx+1,:))/2;
vf(2:2:Mxf,2:2:Myf) = (v(1:Mx,1:My)+v(2:Mx+1,1:My)...
                      + v(1:Mx,2:My+1) + v(2:Mx+1,2:My+1))/4;
end
