%% System variables
syms i_A
syms omega
syms phi
syms x
syms x_dot
syms t

%% Inputs
syms u_A
state = sym('state', [6 1]);

%% Parameters
R_A = 1;
L_A = 1;

k_w = 1;

z = 1;
omega_x = 1;
a = 1;
I_Rinf = 1;

k_m = 1;
eta_G = 1;
J = 1;

r = 1;
G = 1;

b = 1;

c = 1;
d = 1;
Frmin = 1;
Frmax = 2*Frmin;
xmax = 1;

m = 1;
g = 1;
%% algebraic equations
F_c = c * (phi * r/G - x);
F_d = d * (omega * r/G - x_dot);
F_r = sign(x_dot) * (Frmin + x / xmax * (Frmax - Frmin));
F_g = m*g;

M_A = k_m * i_A * eta_G;
M_c = F_c / (r * G);
M_d = F_d / (r * G);
M_L = b * omega + M_c + M_d;

%% ODEs
f_i_A = (u_A * R_A * i_A - k_w * omega) / L_A;

f_omega = (state(1,:) + state(2,:) + state(4,:) + state(5,:) + state(6,:)) * ((M_A - M_L) / (J * eta_G)) + ...
          state(3,:) * 0;
      
f_phi = (state(1,:) + state(2,:) + state(4,:) + state(5,:) + state(6,:)) * omega + ...
        state(3,:) * 0;
    
f_x =  (state(1,:) + state(2,:) + state(3,:)) * x_dot + ...
       (state(4,:) + state(5,:) + state(6,:)) * 0;

f_x_dot =   (state(1,:) + state(2,:) + state(3,:) + state(5,:) + state(6,:)) * ((F_c + F_d - F_r - F_g) * 1/m) + ...
            (state(4,:)) * 0;
        
i_R = ((i_A * a + 1) * I_Rinf * omega)/sqrt(omega_x^2 + omega^2) * sin(z*phi);

y_out = [1, 0, 0, 0, 0, 1];

F_ode = [f_i_A; f_omega; f_phi; f_x; f_x_dot; i_R];
F_ode = matlabFunction(F_ode, 'Vars', {t, [i_A, omega, phi, x, x_dot], u_A, state});

%% algebraic equations
i = i_A + i_R;
