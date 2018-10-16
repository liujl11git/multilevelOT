function [l] = PrimalFunL1(m, dx)
% A function to calculate the 1,1-norm of m.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

l = sum(sum(sum(abs(m))));
l = l * dx * dx;

end
