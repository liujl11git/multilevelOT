function [l] = PrimalFunL2(m, dx)
% A function to calculate the 2,1-norm of m.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

l = sum(sum(sqrt(m(:,:,1).^2 + m(:,:,2).^2)));
l = l * dx * dx;

end
