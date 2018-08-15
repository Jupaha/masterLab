%%  Plot stability region of RK4 
clear i mu;
[X,Y] = meshgrid(-5:0.01:5,-5:0.01:5);
mu = X+i*Y;
R = (mu.^4)/24 + (mu.^3)/6 + (mu.^2)/2 + mu + 1;
Rhat = abs(R);
R_bool = Rhat <= 1;
contour(X,Y,R_bool,'-m')
hold on

mu_min = h_max * lbd;
mu_r = real(mu_min);
mu_i = imag(mu_min);

plot(mu_r, mu_i, 'or', 'LineWidth',10)

% Plot axis
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
grid


