close all
clear all

global blocked;
blocked = false;
global iced;
iced = false;
global phi_switch;
phi_switch = 0;
%global up;
%up = false;
inistate = 5;

y0 = [0 0 0 0 0];

%% Solver
[fct, youtfct, Jf] = getF_window();
[tout,yout,states] = RK4(fct, [0 10], y0, 1e-4, 1e-4, 1e-7, @fensterReset, youtfct, @u_input, inistate);


%% Plot
fig = figure('Position', [-788 -179 789 984]);
hold on
displaynames = ["i_M", "omega", "phi", "x", "x_dot"];
for it=1:size(yout,2)
    ax(it) = subplot(size(yout,2)+1,1,it);
    hold on

    plot(tout, yout(:,it), 'displayName', displaynames(it));
    legend;
end
    ax(size(yout,2)+1) = subplot(size(yout,2)+1,1,size(yout,2)+1);
    hold on
    plot(tout, states, 'displayName', 'states');
    legend;
    linkaxes(ax,'x')
    %% testing
    
    
%% Model Functions
function [F_ode_cells, y_out, Jf] = getF_window()
%% System variables
% 1:i_M; 2:w; 3:phi; 4:x2; 5:x_dot;
syms t_sym   i_M_sym  w_sym  phi_sym  x2_sym  dx2_sym;
y_sym =     [i_M_sym; w_sym; phi_sym; x2_sym; dx2_sym];
syms u_B_sym;
u_sym = [u_B_sym];

R = 0.875;    % motor resistance (0.75)
L = 1500e-6; % motor inductance
k = 0.021;   % motor constant
J = 1.5e-5;  % moment of inertia
b = 2e-5;    % viscouss friction motor (2)

% mechanical parameters of window regulator
m = 8.26;   % equivalent window mass
c = 45e3;    % spring constant of mech. system
d = 310;    % damping constant of mech. system
eta_G = 0.395; % worm gear efficiency
r = 0.025;   % radius of cable drum (0.025)
G = 73;      % worm gear translation ratio

% parameter for variable friction
H = G * (k^2 * eta_G + b * R) / (k * r); % constant for easier use
I_0 = 1.2;
I_up_min = 6;
I_up_max = 7.5;
I_down_min = 3;
I_deltaminmax = I_up_max - I_up_min;
F_r_min = (I_down_min - I_0) * H;
F_r_max = I_deltaminmax * H + F_r_min;

% misc.
g = 9.81;      % gravity
x_max = 0.51;

% ripple parameter
a = 0.595;                   % current influence factor
I_RS = 0.0726;               % peak value current ripple at fast speed
n_x = 20;                     % corner speed
alpha_start = 0.45 * 2 * pi; % start angle
z = 10;                      % commutator segments

di = (u_B_sym - R * i_M_sym - k * w_sym) / L;
F_c = c * (phi_sym * (r / G) - x2_sym);
F_d = d * (w_sym * (r / G) - dx2_sym);
F_r = sign(dx2_sym) *(F_r_min + x2_sym * (F_r_max - F_r_min) / x_max);
F_m = sign(dx2_sym) * m * g;
M_A = k * eta_G * i_M_sym;
M_d = F_d * (r / G);
M_c = F_c * (r / G);
M_b = b * w_sym;
dw = (M_A - M_b - M_d - M_c) / (J * eta_G);
dphi = w_sym;
dx2 = dx2_sym;
ddx2 = (F_c + F_d - F_r - F_m) / m;

% STATE 1; 2;
F_ode = [di; dw; dphi; dx2; ddx2];
F_ode = matlabFunction(F_ode, 'Vars', {t_sym, y_sym, u_sym});

% STATE 3:
F_ode_3 = [di; 0; dphi; dx2; ddx2];
F_ode_3 = matlabFunction(F_ode_3, 'Vars', {t_sym, y_sym, u_sym});

% STATE 4:
F_ode_4 = [di; dw; dphi; 0; 0];
F_ode_4 = matlabFunction(F_ode_4, 'Vars', {t_sym, y_sym, u_sym});

% STATE 5:
dw_5 = (M_A - M_b) / (J * eta_G);
F_ode_5 = [di; dw_5; dphi; 0; 0];
F_ode_5 = matlabFunction(F_ode_5, 'Vars', {t_sym, y_sym, u_sym});

% STATE 6:
F_ode_6 = [di; dw; dphi; 0; 0];
F_ode_6 = matlabFunction(F_ode_6, 'Vars', {t_sym, y_sym, u_sym});

%% ODE
F_ode_cells = {F_ode, F_ode, F_ode_3, F_ode_4, F_ode_5, F_ode_6};

n = @(omega)omega ./ (2 .* pi);
y_out = @(y_state)[   y_state(:,1) + (((a .* abs(y_state(:,1)) + 1) .* I_RS .* n(y_state(:,2))) ./ sqrt(n(y_state(:,2)) .* n(y_state(:,2)) + n_x.^2) .* sin((y_state(:,3) + alpha_start) .* z)),...
            y_state(:,2),...
            y_state(:,3),...
            y_state(:,4),...
            y_state(:,5)];

Jf = cell(1,length(F_ode_cells));
for i=1:5
    fct = F_ode_cells{i};
    Jf_sym = jacobian(fct(t_sym,y_sym,u_sym),y_sym);
    Jf{i} = matlabFunction(Jf_sym,'Vars',{t_sym,y_sym,u_sym});
end
end
function [y_out, state] = fensterReset(t,y_in,u,state)
i_M     = y_in(1);
w       = y_in(2);
phi     = y_in(3);
x2       = y_in(4);
dx2      = y_in(5);
u_B     = u(1);


R = 0.875;    % motor resistance (0.75)
L = 1500e-6; % motor inductance
k = 0.021;   % motor constant
J = 1.5e-5;  % moment of inertia
b = 2e-5;    % viscouss friction motor (2)

% mechanical parameters of window regulator
m = 8.26;   % equivalent window mass
c = 45e3;    % spring constant of mech. system
d = 310;    % damping constant of mech. system
eta_G = 0.395; % worm gear efficiency
r = 0.025;   % radius of cable drum (0.025)
G = 73;      % worm gear translation ratio

% parameter for variable friction
H = G * (k^2 * eta_G + b * R) / (k * r); % constant for easier use
I_0 = 1.2;
I_up_min = 6;
I_up_max = 7.5;
I_down_min = 3;
I_down_max = 4.5;
I_deltaminmax = I_up_max - I_up_min;
I_m = I_down_min - I_0;
I_delta_updown = I_up_min - I_down_min;
I_n = (I_up_min + I_up_max + I_down_min + I_down_max) / 4;
F_r_min = (I_down_min - I_0) * H;
F_r_max = I_deltaminmax * H + F_r_min;

% misc.
g = 9.81;      % gravity
x_max = 0.51;
global iced;
global phi_switch;
global blocked;
%global up;
phi_loose = 61;


% ripple parameter
a = 0.595;                   % current influence factor
I_RS = 0.0726;               % peak value current ripple at fast speed
n_x = 20;                     % corner speed
alpha_start = 0.45 * 2 * pi; % start angle
z = 10;                      % commutator segments

F_break = H * 0.78 * u_B / R;
h_B = 0.1;

switch state 
    case 1
        F_c = c * (phi * (r / G) - x2);
        F_d = d * (w * (r / G) - dx2);
        M_A = k * eta_G * i_M;
        M_d = F_d * (r / G);
        M_c = F_c * (r / G);
        M_b = b * w;
        M_L = M_d + M_c + M_b;

        if abs(M_A) < abs(M_L) && sign(M_L) ~= sign(w)
            state = 3;
        elseif x2 > x_max
            state = 4;
            x2 = x_max;
            %up = true;
        elseif blocked && x2 > x_max - h_B
            state = 4;
            x2 = x_max - h_B;
            %up = true;
        end

    case 2
        F_c = c * (phi * (r / G) - x2);
        F_d = d * (w * (r / G) - dx2);
        M_A = k * eta_G * i_M;
        M_d = F_d * (r / G);
        M_c = F_c * (r / G);
        M_b = b * w;
        M_L = M_d + M_c + M_b;
        
        if x2 < 0
            state = 4;
            x2 = 0;
            %up = false;
        elseif abs(M_A) < abs(M_L) && sign(M_A) ~= sign(w)
            state = 3;
        end
        
    case 3
        F_c = c * (phi * (r / G) - x2);
        F_d = d * (w * (r / G) - dx2);
        M_A = k * eta_G * i_M;
        M_d = F_d * (r / G);
        M_c = F_c * (r / G);
        M_b = b * w;
        M_L = M_d + M_c + M_b;
        
        n = w ./ (2 .* pi);
        i = i_M + (((a .* abs(i_M) + 1) .* I_RS .* n ./ sqrt(n .* n + n_x.^2)) .* sin((phi + alpha_start) .* z));
       
        
        if abs(M_A) > abs(M_L) && sign(M_A) == sign(M_L) && i > 0
            state = 1;
        elseif sign(M_A) ~= sign(M_L) && M_A > 0
            state = 1;
        elseif abs(M_A) > abs(M_L) && sign(M_A) == sign(M_L) && i < 0
            state = 2;
        elseif abs(M_A) > abs(M_L) && sign(M_A) ~= sign(M_L) && i < 0
            state = 2;
        else
            w = 0;
            if x2 > x_max
                x2 = x_max;
                %up = true;
            elseif x2 < 0
                    x2 = 0;
                    %up = false;
            end
        end
    case 4
        F_c = c * (phi * (r / G) - x2);
        F_d = d * (w * (r / G) - dx2);
        M_A = k * eta_G * i_M;
        M_d = F_d * (r / G);
        M_c = F_c * (r / G);
        M_b = b * w;
        M_L = M_d + M_c + M_b;
        
        n = w ./ (2 .* pi);
        i = i_M + (((a .* abs(i_M) + 1) .* I_RS .* n ./ sqrt(n .* n + n_x.^2)) .* sin((phi + alpha_start) .* z));
        
        
        if abs(u_B) < 1 && abs(i_M) < 0.1
            if x2 <= 0
                phi = -phi_loose;
                phi_switch = 0;
                state = 5;
                %up = false;
            elseif x2 >= x_max
                phi = (x_max * (G / r)) + phi_loose;
                phi_switch = (x_max * (G / r));
                state = 5;
                %up = false;
            elseif blocked && x2 >= x_max - h_B
                phi = ((x_max-h_B) * (G / r)) + phi_loose;
                phi_switch = ((x_max-h_B) * (G / r));
                state = 5;
                %up = false;
            end
        elseif abs(M_A) < abs(M_L) && sign(M_A) == sign(M_L) && M_L ~= 0 && u == 0 && abs(i) < 0.1
            state = 3;
            %up = false;
        elseif x2 > 0 && x2 <= x_max && dx2 > 0
            %elseif ~up &&  dx2 > 0
            state = 1;
            %up = false;
        elseif x2 > 0 && x2 <= x_max && dx2 < 0
            %elseif up && dx2 < 0
            state = 2;
            %up = false;
        elseif sign(M_A) ~= sign(M_L) && M_A > 0
            state = 1;
            %up = false;
        else
            dx2 = 0;
        end
                

    case 5 
        if w > 0 && phi > phi_switch
            state = 1;
        end
        if w < 0 && phi < phi_switch
            if iced
                state = 6;
            else
                state = 2;
            end
        end
    case 6
        F_c = c * (phi * (r / G) - x2);
        F_d = d * (w * (r / G) - dx2);
        F_r = sign(dx2_sym) *(F_r_min + x2_sym * (F_r_max - F_r_min) / x_max);    
        
        if F_c + F_d - F_r < -F_break
            state = 2;
            iced = False;
        else
            dx2 = 0;
        end
    otherwise
        n = w ./ (2 .* pi);
        i = i_M + (((a .* abs(i_M) + 1) .* I_RS .* n ./ sqrt(n .* n + n_x.^2)) .* sin((phi + alpha_start) .* z));
        
        if phi > 0 && i > 0
            state = 1;
        elseif phi < phi_max && i < 0
            state = 2;
        elseif (phi < 0 && i > 0) || (phi > phi_max && i < 0)
            state = 5;
        end
end
        
y_out(1) = i_M;
y_out(2) = w;
y_out(3) = phi;
y_out(4) = x2;
y_out(5) = dx2;    
end
function u = u_input(t)
    if (t < 6) 
        u(1,:) = 12.5;
    else
        u(1,:) = 0;
    end
end

%% Solver Functions
function [t,y_out,states] = RK4(f, dT, y0, h_start, h_min, tau, reset_y, get_y_out, get_u, initState)
% Klassisches Runge-Kutta verfahren (m=4)
disp('RK4')

tic
t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);
states(t_count) = initState;

% RK coefficients
a = [ 0,  0, 0, 0;...
    1/2,  0, 0, 0;...
      0,1/2, 0, 0;...
      0,  0, 1, 0];
b = [1/6, 1/3, 1/3, 1/6];
c = [0; 1/2; 1/2; 1];

while t(t_count)<dT(2)
    u_in = get_u(t(t_count));
    y(t_count+1,:) = RK_step(a,b,c,f{states(t_count)},t(t_count),y(t_count,:),u_in,h); 
    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
    
    [y(t_count,:), states(t_count)] = reset_y(t(t_count),y(t_count,:),u_in,states(t_count-1));
        
end

y_out = get_y_out(y);

el_ode = toc;
disp(['-Steps:           ' num2str(length(t))])
disp(['-Time:            ' num2str(el_ode)])
disp(' ')
t=t.';
end
function [t,y_out,states] = RK_imp6(f, Jf, dT, y0, h_start, h_min, tau, reset_y, get_y_out, get_u, initState)
% IMPLICIT RUNGE-KUTTA ALGORITHM USING NEWTON-RAPHSON METHOD by Andr´es L. Granados M
disp('RK_imp6')
tic

n = 1;
t(n) = dT(1);
h = max(h_start, h_min);
y(n,:) = y0(:);
states(n) = initState;

% explicit RK coefficients
a_exp = [(5-sqrt(15))/10,                0,         0;
        -(3+sqrt(15))/4,    (5+sqrt(15))/4,         0;
        3*(4+sqrt(15))/5,    -(35+9*sqrt(15))/10,    2*(4+sqrt(15))/5];

% implicit RK coefficients
a = [5/36               (10-3*sqrt(15))/45  (25-6*sqrt(15))/180;
    (10+3*sqrt(15))/72  2/9                 (10-3*sqrt(15))/72
    (25+6*sqrt(15))/180 (10+3*sqrt(15))/45  5/36];

c = [(5-sqrt(15))/10 1/2 (5+sqrt(15))/10];

b = [5/18 4/9 5/18];

    
omega = 1.0;

N = length(y0);

I = eye(N);

Jg = @(Jf,tn,yn,un,hn,knsum)[   hn*a(1,1)*(Jf(tn + c(1) * hn, (yn + hn * knsum(1,:)).', un))-I,         hn*a(1,2)*(Jf(tn + c(2) * hn, (yn + hn * knsum(1,:)).', un)),           hn*a(1,3)*(Jf(tn + c(3) * hn, (yn + hn * knsum(1,:)).', un));
                                hn*a(2,1)*(Jf(tn + c(1) * hn, (yn + hn * knsum(2,:)).', un)),           hn*a(2,2)*(Jf(tn + c(2) * hn, (yn + hn * knsum(2,:)).', un))-I,         hn*a(2,3)*(Jf(tn + c(3) * hn, (yn + hn * knsum(2,:)).', un));
                                hn*a(3,1)*(Jf(tn + c(1) * hn, (yn + hn * knsum(3,:)).', un)),           hn*a(3,2)*(Jf(tn + c(2) * hn, (yn + hn * knsum(3,:)).', un)),           hn*a(3,3)*(Jf(tn + c(3) * hn, (yn + hn * knsum(3,:)).', un))-I];

                   
while t(n)<dT(2)
    f_state = f{states(n)};
    Jf_state = Jf{states(n)};
    u_in = get_u(t(n));
    
    % estimation of initial k
    k(1,:) = f_state(t(n), y(n,:).', u_in).';
    k(2,:) = f_state(t(n) + c(2) * h, (y(n,:) + h * a_exp(2,1) * k(1,:)).', u_in).';    
    k(3,:) = f_state(t(n) + c(3) * h, (y(n,:) + h * (a_exp(3,1) * k(1,:)) + a_exp(3,2) * k(2,:)).', u_in).'; 
    
    err = 1;
    while max(err)>1e-2
    ksum = [(a(1,1) * k(1,:) + a(1,2) * k(2,:) + a(1,3) * k(3,:));
            (a(2,1) * k(1,:) + a(2,2) * k(2,:) + a(2,3) * k(3,:));
            (a(3,1) * k(1,:) + a(3,2) * k(2,:) + a(3,3) * k(3,:))];
    
    g = [f_state(t(n) + c(1) * h, (y(n,:) + h * ksum(1,:)).', u_in)-k(1,:).';
         f_state(t(n) + c(2) * h, (y(n,:) + h * ksum(2,:)).', u_in)-k(2,:).';
         f_state(t(n) + c(3) * h, (y(n,:) + h * ksum(3,:)).', u_in)-k(3,:).'];
      
    Jgn = Jg(Jf_state,t(n),y(n,:),u_in,h,ksum);
    
    dk_srt = linsolve(Jgn,-1*g);
    dk = [dk_srt(1:N).';
          dk_srt(N+1:2*N).';
          dk_srt(2*N+1:3*N).'];
    
    k = k + omega * dk;
    
    err = abs(omega*h*[c(1)*norm(dk(1,:)); c(2)*norm(dk(2,:)); c(3)*norm(dk(3,:))]);
    end
    
    y(n+1,:) = y(n,:) + h * (b(1)*k(1,:) + b(2)*k(2,:) + b(3)*k(3,:));
    n = n+1;
    t(n) = t(n-1) + h;
    
    [y(n,:), states(n)] = reset_y(t(n),y(n,:),u_in,states(n-1));
end

y_out = get_y_out(y);

el_ode = toc;
disp(['-Steps:           ' num2str(length(t))])
disp(['-Time:            ' num2str(el_ode)])
disp(' ')

t=t.';
end
function [yn1] = RK_step(a,b,c,f,t,yn,u,h)
    s = length(c);
    k(1:s,1:length(yn)) = 0;
    for i=1:s
        k(i,:)   = f(t + c(i) * h, (yn + h * (a(i,:)*k)).', u).';
    end
    yn1 = yn + h * b * k;    
end 