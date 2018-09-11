close all
%% System variables
syms i_A
syms omega
syms phi
syms x
syms x_dot
syms t

%% Inputs
syms u_A;

%% Parameters
R_A = 0.75;
L_A = 1.5e-3;
k_m = 0.021; k_w = k_m;
J = 15e-6;
a = 20e-6;
I_Rinf = 0.0726;
omega_x = 125;
z = 10;
G = 73;
eta_G = 0.395;
xmax = 0.49;
m = 8.26;
r = 0.025;
c = 45e3;
d = 310;
Phi_L = 55;
vereist = false;
F_b = 350;
blockiert = false;
h_b = 0.1;
Frmin = 47;
Frmax = 87;

b = 1e-5;
g = 9.81;

phi_s = 0;

%% algebraic equations
F_c = c * (phi * r/G - x);
F_d = d * (omega * r/G - x_dot);
F_r = sign(x_dot) * (Frmin + abs(x / xmax) * (Frmax - Frmin));
F_g = m*g;

M_A = k_m * i_A * eta_G;
M_c = F_c / (r * G);
M_d = F_d / (r * G);
M_L = -b * omega - M_c - M_d;

u_i = k_w * omega;

%% ODEs
f_i_A = (u_A - R_A * i_A - u_i) / L_A;
f_omega = ((M_A + M_L) / (J * eta_G)); %((M_A - b*omega) / (J * eta_G))
f_phi = omega;  
f_x =  x_dot;
f_x_dot =   ((F_c + F_d - F_r - F_g) * 1/m); %-g    

F_ode = [f_i_A; f_omega; f_phi; f_x; f_x_dot];
F_ode = matlabFunction(F_ode, 'Vars', {t, [i_A; omega; phi; x; x_dot], u_A});

i_R = ((i_A * a + 1) * I_Rinf * omega)/sqrt(omega_x^2 + omega^2) * sin(z*phi);
i = i_A + i_R;
%% Solver
y0 = [0; 0; 0; 0; 0];
[tout,yout] = RK4(F_ode, [0 2], y0, 1e-4, 1e-4, 1e-7, @fensterReset);

%% Plot
fig = figure('Position', [-788 -179 789 984]);
hold on
displaynames = ["i_A", "omega", "phi", "x", "x_dot"];
for it=1:size(yout,2)
    subplot(size(yout,2),1,it)
    hold on

    plot(tout, yout(:,it), 'displayName', displaynames(it));
    legend;
end
%% Plot
% subplot(2,1,1)
% plot(tout, yout(:,1))
% subplot(2,1,2)
% plot(tout, yout(:,2))

function [y_out, state] = fensterReset(t, y_in, state)
i_A     = y_in(1);
omega   = y_in(2);
phi     = y_in(3);
x       = y_in(4);
x_dot   = y_in(5);

if false
% R_A = 0.75;
% L_A = 1.5e-3;
% k_m = 0.021; k_w = k_m;
% J = 15e-6;
% a = 20e-6;
% I_Rinf = 0.0726;
% omega_x = 125;
% z = 10;
% G = 73;
% eta_G = 0.395;
% xmax = 0.49;
% m = 8.26;
% r = 0.025;
% c = 45e3;
% d = 310;
% Phi_L = 55;
% vereist = false;
% F_b = 350;
% blockiert = false;
% h_b = 0.1;
% Frmin = 47;
% Frmax = 87;
% 
% b = 1e-5;
% g = 9.81;
% 
% phi_s = 0;
% 
% x_anschlag = 0;
% 
% F_c = c * (phi * r/G - x);
% F_d = d * (omega * r/G - x_dot);
% F_r = sign(x_dot) * (Frmin + abs(x / xmax) * (Frmax - Frmin));
% F_g = m*g;
% 
% M_A = k_m * i_A * eta_G;
% M_c = F_c / (r * G);
% M_d = F_d / (r * G);
% M_L = -b * omega - M_c - M_d;
% 
% u_i = k_w * omega;
% 
% cond13 = abs(M_A) < abs(M_L) && sign(M_L) ~= sign(omega);
% cond14 = x >= xmax;
% cond23 = cond13;
% cond24 = x <= 0;
% cond31 = (abs(M_A) > abs(M_L) && sign(M_A) == sign(M_L) && i_A > 0) || ...
%          (sign(M_A) ~= sign(M_L) && M_A > 1e-4);
% cond32 = (abs(M_A) > abs(M_L) && sign(M_A) == sign(M_L) && i_A < 0) || ...
%          (sign(M_A) ~= sign(M_L) && M_A < -1e-4);
% cond41 = sign(M_A) ~= sign(M_L) && M_A > 1e-4;
% cond42 = sign(M_A) ~= sign(M_L) && M_A < -1e-4;
% cond43 = blockiert;
% cond45 = false;
% cond51 = omega > 0 && phi > phi_s;
% cond52 = omega < 0 && phi < phi_s && ~vereist;
% cond56 = omega < 0 && phi < phi_s && vereist;
% cond62 = F_c + F_d - F_r < -F_b;
% 
% switch state
%     case 1
%         if cond13
%             state = 3;
%         elseif cond14
%             state = 4;
%             % Fenster am oberen Anschlag
%             x_anschlag = xmax;
%         end
%     case 2
%         if cond23
%             state = 3;
%         elseif cond24
%             state = 4;
%             % Fenster am unteren Anschlag
%             x_anschlag = 0;
%         end
%     case 3
%         if cond31
%             state = 1;
%         elseif cond32
%             state = 2;
%         end
%     case 4
%         if cond41
%             state = 1;
%         elseif cond42
%             state = 2;
%         elseif cond43
%             state = 3;
%         elseif cond45
%             state = 5;
%         end
%     case 5
%         if cond51
%             state = 1;
%         elseif cond52
%             state = 2;
%         elseif cond56
%             state = 6;
%         end
%     case 6
%         if cond62
%             state = 2;
%         end
%     otherwise
%         state = 4;
% end
% 
% switch state
%     case 1
%     case 2
%         vereist = 0;
%     case 3
%         omega = 0;
%     case 4
%         x = x_anschlag;
%         x_dot = 0;
%     case 5
% %         F_c = 0;
% %         F_d = 0;
% %         F_r = 0;
% %         x_dot = 0;
%     case 6
%         x_dot = 0;
% end
end
y_out(1) = i_A;
y_out(2) = omega;
y_out(3) = phi;
y_out(4) = x;
y_out(5) = x_dot;

end
function [t,y] = RK4(f, dT, y0, h_start, h_min, tau, resetFunction)
% Klassisches Runge-Kutta verfahren (m=4)
disp('RK4')
tic

t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);
state = 0;

% RK coefficients
a = [ 0,  0, 0, 0;...
    1/2,  0, 0, 0;...
      0,1/2, 0, 0;...
      0,  0, 1, 0];
b = [1/6, 1/3, 1/3, 1/6];
c = [0; 1/2; 1/2; 1];


while t(t_count)<dT(2)         
    y(t_count+1,:) = RK_step(a,b,c,f,y(t_count,:),t(t_count),h);
    
    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
    
    if exist('resetFunction', 'var')
      [y(t_count,:) state] = resetFunction(t(t_count), y(t_count,:), state);
    end
end
el_ode = toc;
disp(['-Steps:           ' num2str(length(t))])
disp(['-Time:            ' num2str(el_ode)])
disp(' ')
t=t.';
end
function [yn1] = RK_step(a,b,c,f,yn,t,h)
    s = length(c);
    k(1:s,1:length(yn)) = 0;
    for i=1:s
        k(i,:)   = f(t + c(i) * h, (yn + h * (a(i,:)*k)).', 240).';
    end
    yn1 = yn + h * b * k;    
end 