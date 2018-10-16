function [phi] = Div(m,dx) 
% A function to calculate the divergence of flux m
% Author: Jialin Liu (liujl11@math.ucla.edu) 
% Author: Wuchen Li (wcli@math.ucla.edu) 
% Modified: 2018-10-10

Mx = size(m,1);
My = size(m,2);
phi = zeros(Mx,My);
phi(2:Mx,:) = phi(2:Mx,:) + (m(2:end,:,1) - m(1:end-1,:,1) )/dx;
phi(1,:) = phi(1,:) + m(1,:,1)/dx;
phi(:,2:My) = phi(:,2:My) + (m(:,2:end,2) - m(:,1:end-1,2))/dx;
phi(:,1) = phi(:,1) + m(:,1,2)/dx;
end
