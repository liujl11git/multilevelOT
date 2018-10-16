function [h,rho0,rho1,x,y] = ReadIMG(img1,img2)
% A function to create a OT problem from two images.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

img1 = double(imread(img1));
img2 = double(imread(img2));

if size(img1,3)>1, img1 = mean(img1,3); end
if size(img2,3)>1, img2 = mean(img2,3); end

[Mx1, My1] = size(img1);
[Mx2, My2] = size(img2);
M = max([Mx1 My1 Mx2 My2]);

Mpower = ceil(log(M-1)/log(2));
M = 2^Mpower + 1;

h=1/(M-1);
x = 0:h:1;
y = 0:h:1;

% The following operation is to add a zero boundary
% The reason why we do this is:
% In the grid system, if there are 256*256 "squares",
% there are 257*257 nodes.
rho0 = padarray(img1, [M-Mx1, M-My1], 0, 'post');
rho1 = padarray(img2, [M-Mx2, M-My2], 0, 'post');

min1 = min(rho0(:));
min2 = min(rho1(:));
if min1 < 0, rho0 = rho0 - min1; end % must be non-negative
if min2 < 0, rho1 = rho1 - min2; end

rho0 = rho0 / sum(rho0(:));
rho1 = rho1 / sum(rho1(:));
rho0 = rho0 / h / h;
rho1 = rho1 / h / h;
end
