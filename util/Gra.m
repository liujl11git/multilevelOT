function [m] = Gra(phi,dx)
% A function to calculate the gradient of potential phi.
% Author: Jialin Liu (liujl11@math.ucla.edu) 
% Author: Wuchen Li (wcli@math.ucla.edu) 
% Modified: 2018-10-10

Mx = size(phi,1);
My = size(phi,2);
m = zeros(Mx,My,2);
m(1:Mx-1,:,1) = (phi(2:Mx,:) - phi(1:Mx-1,:))/dx;
m(:,1:My-1,2) = (phi(:,2:My) - phi(:,1:My-1))/dx;
end
