%%  Plot stability region of RK4 
clearvars -except eigV h_eigmax
clear i;
[X,Y] = meshgrid(-10:0.005:10,-10:0.005:10);
m = X+i*Y;

% stability function R see: Stability of Runge-Kutta Methods
R = (m.^4)/24 + (m.^3)/6 + (m.^2)/2 + m + 1;
%R = m.^7./24000 + (161.*m.^6)./120000 + (195398609418584047.*m.^5)./21374506043965440000 + (118747255799808017.*m.^4)./2849934139195392000 + (419573637159321617.*m.^3)./2517441822955929600 + (559431516212428817.*m.^2)./1118863032424857600 + (2983634753132953651.*m)./2983634753132953600 + 1;

Rhat = abs(R);
R_bool = Rhat <= 1;
%%
contour(X,Y,R_bool,'-m')
hold on

plot(real(eigV), imag(eigV), '*')
%%
% x = -5;
% y = 7;
% mu_min = x+y*i;
% plot(x, y, 'or', 'LineWidth',10)
% line([x 0], [y 0]);  %line
% F_bool = abs(X+i*Y)<=abs(mu_min);
% contour(X,Y,F_bool,'-r')
% a_sqr = (abs(x)-abs(real(scaled_mu)))^2;
% b_sqr = (abs(y)-abs(imag(scaled_mu)))^2;
% r = sqrt((a_sqr)+(b_sqr));
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% h = plot(xunit, yunit);
%%
% Plot axis
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
grid


