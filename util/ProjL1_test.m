% A script to test "projL1.m"
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

D = 0.05 * randn(10,10,2);

nnz(abs(D(:,:,1))+abs(D(:,:,2))<=0.3)/100

D1 = ProjL1(D,1);

nnz(abs(D1(:,:,1))+abs(D1(:,:,2))<=0.3)/100

dotsx = reshape(D(:,:,1),[100 1]);
dotsy = reshape(D(:,:,2),[100 1]);

dotsx1 = reshape(D1(:,:,1),[100 1]);
dotsy1 = reshape(D1(:,:,2),[100 1]);

samplep = 1:100;

polyx = [0 1 0 -1];
polyy = [1 0 -1 0];
fill(polyx,polyy,'y');
hold on;
plot(dotsx(samplep),dotsy(samplep),'b*');
plot(dotsx1(samplep),dotsy1(samplep),'ro');
hold off;
