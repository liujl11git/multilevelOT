function [l] = PrimalFunLinf(m, dx)
% A function to calculate the inf,1-norm of m.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

l = sum(sum( max(abs(m(:,:,1)), abs(m(:,:,2))) ));
l = l * dx * dx;

end
