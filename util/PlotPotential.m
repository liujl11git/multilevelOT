function [] = PlotPotential(x,y,rho)
% A function to plot the potential rho.
% Author: Jialin Liu (liujl11@math.ucla.edu) Modified: 2018-10-10

hx = x(2) - x(1);
hy = y(2) - y(1);

figure ('Units', 'pixels', 'Position', [0 0 400 400]) ;
imagesc('XData',x,'YData',y,'CData',(rho)','CDataMapping','scaled');

set(gca, 'YLim', [-hy/2 1+hy/2], 'XLim',[-hx/2 1+hx/2]);
camroll(-90);

title('The dual solution Phi');
end
