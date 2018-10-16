function [] = PlotFlow(x,y,m,rho)
% A function to plot the flow m with potential rho
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

hx = x(2) - x(1);
hy = y(2) - y(1);

figure ('Units', 'pixels', 'Position', [0 0 400 400]) ;
imagesc('XData',x,'YData',y,'CData',(rho)','CDataMapping','scaled');
hold on;

mx = m(:,:,1);
mx(2:end,:) = (m(1:end-1,:,1)+m(2:end,:,1))/2;

my = m(:,:,2);
my(:,2:end) = (m(:,1:end-1,2)+m(:,2:end,2))/2;

if size(m,1)>20
[x1,y1] = meshgrid(0:(1/32):1,0:(1/32):1);
msmooth = zeros(size(m));
msmooth(:,:,1) = mx;
msmooth(:,:,2) = my;
while size(msmooth,1)>33
    msmooth = restrict_m(msmooth);
end

quiver(x1,y1, (msmooth(:,:,1))', (msmooth(:,:,2))','k');
hold off;
end

if size(m,1)<=20
[x1,y1] = meshgrid(x,y);
msmooth = zeros(size(m));
msmooth(:,:,1) = mx;
msmooth(:,:,2) = my;
quiver(x1,y1, (msmooth(:,:,1))', (msmooth(:,:,2))','k');
hold off;
end

set(gca, 'YLim', [-hy/2 1+hy/2], 'XLim',[-hx/2 1+hx/2]);
camroll(-90);

title('Optimal FLux m');
end
