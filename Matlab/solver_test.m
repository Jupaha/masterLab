%% Numerische Lösung durch verschiedene Solver
    close all
    clear all
    disp(' ')
    disp('%%%')

%% Set Up ODE function
    % Case 1: Dahlquist
    % Case 2: Stiff 2x2-System
    % Case 3: halogen ODE
    % Case 4: Orbiting
    ode_case = 3;
    lambda = -3;                    % Lambda of ODE
    dT = [0 10];                    % Simulation time [t_start t_end]
    tol = 1e-7;                     % Tolerance       
    [fct, fct_sol, y0] = getfct(ode_case,lambda);
    

%% Analyze ODE
% Range of variables
    r = [0, 100; 300 10000];
% define steps in ranges
    n = 10;    

% Compute all eigenvalues in range
eigV_all = comp_eigenvalues(fct,r,n);
eigV = unique(eigV_all);

lmbd_min = min(eigV);

a_RK4 = [0  0  0 0;
    1/2 0  0 0;
    0  1/2 0 0;
    0   0  1 0];
b_RK4 = [1/6 1/3 1/3 1/6];
c_RK4 = [0 1/2 172 1];

h_max_RK4 = get_h_max(a_RK4,b_RK4,c_RK4,lmbd_min);
    
 
%% Time
    h_start = 0.002;
    h_min = h_start;
    
%% Solver selection:

% RK_ODE45, RK_ODE4, RK_imp6, RK_imp6_auto, RK_xorder, RK_DP, RK4_auto,
% RK4, FWE_auto, FWE2, FWE, BWE_fp
solver_selection = ["RK_ODE45", "RK_ODE4", "RK_imp6", "RK_imp6_auto", "RK_xorder", "RK_DP", "RK4_auto", "RK4", "FWE_auto", "FWE2", "FWE", "BWE_fp"];
out = runSolver(solver_selection, fct, dT, y0, h_start, h_min, tol);
    
%% Plot numerical solutions
plotfig = plotNumeric(solver_selection, ode_case, out, length(y0));

%% Plot Exact Solution
if ~isreal(fct_sol)
figure(plotfig)
pT = dT(1):1/((dT(2)-dT(1))*1000):dT(2);
solplot = fct_sol(pT).';
for i=1:length(y0)
subplot(length(y0),1,i)
hold on
plot(pT, solplot(:,i), '--r', 'DisplayName', 'Exact Solution');
legend('Location', 'northeast')
end
end


%% Absolute Error
%     errfig = figure('Position', [0 0 770 800]);
%     hold on
% if exist('fct_sol')
%     abserr_ode = abserr(fct_sol, t_ode, y_ode);
%     %abserr_RKni = abserr(fct_sol, t_RKni, y_RKni);
%     abserr_RKnia = abserr(fct_sol, t_RKnia, y_RKnia);
%     %abserr_FWE = abserr(fct_sol, t_fwe, y_fwe);
% %     abserr_RK4a2 = abserr(fct_sol, t_RK4a2, y_RK4a2);
%     abserr_RKDP = abserr(fct_sol, t_RKDP, y_RKDP);
% 
%     plot(t_ode, abserr_ode, 'DisplayName', 'AbsERR ODE45')
%     plot(t_RKnia, abserr_RKnia, 'DisplayName', 'AbsERR RKnia ')
%     plot(t_RKDP, abserr_RKDP, 'DisplayName', 'AbsERR RKDP')
% end
% 
% legend
%% Error Estimation
% figure(errfig)
% hold on
% plot(t_RK4a2,errest_RK4a2, '--*', 'DisplayName', 'Fehlerschätzung RK4')
% plot(t_RKDP,errest_RKDP, '--*', 'DisplayName', 'Fehlerschätzung RKDP')
%% Step size
% stepfig = figure('Position', [670 0 770 800]);
% t_out = t_RKDP;
% semilogy(t_out(1:end-1), diff(t_out), '-x', 'DisplayName', ['steps: ' num2str(length((t_out)))])
% %axis([0 1 1 1e-8])
% legend
%% SOLVER
function [t,y_out] = RK_ODE45(f, dT, y0, h_start, h_min, tau)
ode_opts = odeset('RelTol',tau, 'AbsTol', tau);
tic
[t,y_out] = ode45(f, dT, y0, ode_opts);
el_ode = toc;
disp(['Steps of RK_ODE45:           ' num2str(length(t))])
disp(['Time of RK_ODE45:            ' num2str(el_ode)])
disp(' ')
end
function [t,y_out] = RK_ODE4(f, dT, y0, h_start, h_min, tau)
h = max(h_start, h_min);
tic
t = dT(1):h:dT(2);
y_out = ode4(f,t,y0);
el_ode = toc;
disp(['Steps of RK_ODE4:            ' num2str(length(t))])
disp(['Time of RK_ODE4:             ' num2str(el_ode)])
disp(' ')
end
function [t,y_out] = RK_imp6(f, dT, y0, h_start, h_min, tau)
% IMPLICIT RUNGE-KUTTA ALGORITHM USING NEWTON-RAPHSON METHOD by Andr´es L. Granados M

syms t_sym;
y_sym = sym('y_sym', [1,length(y0)]);
Jf_sym = jacobian(f(t_sym,y_sym));
Jf = matlabFunction(Jf_sym,'Vars',{t_sym,y_sym});

tic

n = 1;
t(n) = dT(1);
h = h_start;
y_out(n,:) = y0(:);

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

Jg = @(tn,yn,hn,knsum)[ hn*a(1,1)*(Jf(tn + c(1) * hn, yn + hn * knsum(1,:)))-I,    hn*a(1,2)*(Jf(tn + c(2) * hn, yn + hn * knsum(1,:))),           hn*a(1,3)*(Jf(tn + c(3) * hn, yn + hn * knsum(1,:)));
                        hn*a(2,1)*(Jf(tn + c(1) * hn, yn + hn * knsum(2,:))),           hn*a(2,2)*(Jf(tn + c(2) * hn, yn + hn * knsum(2,:)))-I,    hn*a(2,3)*(Jf(tn + c(3) * hn, yn + hn * knsum(2,:)));
                        hn*a(3,1)*(Jf(tn + c(1) * hn, yn + hn * knsum(3,:))),           hn*a(3,2)*(Jf(tn + c(2) * hn, yn + hn * knsum(3,:))),           hn*a(3,3)*(Jf(tn + c(3) * hn, yn + hn * knsum(3,:)))-I];

                   
while t(n)<dT(2)
    
    % estimation of initial k
    k(1,:) = f(t(n), y_out(n,:)).';
    k(2,:) = f(t(n) + c(2) * h, y_out(n,:) + h * a_exp(2,1) * k(1,:)).';    
    k(3,:) = f(t(n) + c(3) * h, y_out(n,:) + h * (a_exp(3,1) * k(1,:) + a_exp(3,2) * k(2,:))).'; 
    
    err = 1;
    while max(err)>1e-12
    ksum = [(a(1,1) * k(1,:) + a(1,2) * k(2,:) + a(1,3) * k(3,:));
            (a(2,1) * k(1,:) + a(2,2) * k(2,:) + a(2,3) * k(3,:));
            (a(3,1) * k(1,:) + a(3,2) * k(2,:) + a(3,3) * k(3,:))];
    
    g = [f(t(n) + c(1) * h, y_out(n,:) + h * ksum(1,:))-k(1,:).';
         f(t(n) + c(2) * h, y_out(n,:) + h * ksum(2,:))-k(2,:).';
         f(t(n) + c(3) * h, y_out(n,:) + h * ksum(3,:))-k(3,:).'];
      
    Jgn = Jg(t(n),y_out(n,:),h,ksum);
    
    dk_srt = linsolve(Jgn,-1*g);
    dk = [dk_srt(1:N).';
          dk_srt(N+1:2*N).';
          dk_srt(2*N+1:3*N).'];
    
    k = k + omega * dk;
    
    err = abs(omega*h*[c(1)*norm(dk(1,:)); c(2)*norm(dk(2,:)); c(3)*norm(dk(3,:))]);
    end
    
    y_out(n+1,:) = y_out(n,:) + h * (b(1)*k(1,:) + b(2)*k(2,:) + b(3)*k(3,:));
    n = n+1;
    t(n) = t(n-1) + h;
end
el_RKni = toc;
disp(['Steps of RK_imp6:            ' num2str(length(t))])
disp(['Time of RK_imp6:             ' num2str(el_RKni)])
%disp(['No. rejects:                ' num2str(rejects)])
disp(' ')

t=t.';
end
function [t,y_out] = RK_imp6_auto(f, dT, y0, h_start, h_min, tau)
% IMPLICIT RUNGE-KUTTA ALGORITHM USING NEWTON-RAPHSON METHOD by Andr´es L. Granados M
% Erweitert mit adaptive step size control

syms t_sym;
y_sym = sym('y_sym', [1,length(y0)]);
Jf_sym = jacobian(f(t_sym,y_sym));
Jf = matlabFunction(Jf_sym,'Vars',{t_sym,y_sym});

tic

rejects = 0;
mu = 2;
rho = 0.8;
n = 1;
t(n) = dT(1);
h = h_start;
y_out(n,:) = y0(:);

% explicit RK coefficients
a_exp = [(5-sqrt(15))/10,                0,         0;
        -(3+sqrt(15))/4,    (5+sqrt(15))/4,         0;
        3*(4+sqrt(15))/5,    -(35+9*sqrt(15))/10,    2*(4+sqrt(15))/5];

p = 6; % RK of 6th order
% implicit RK coefficients
a = [5/36               (10-3*sqrt(15))/45  (25-6*sqrt(15))/180;
    (10+3*sqrt(15))/72  2/9                 (10-3*sqrt(15))/72
    (25+6*sqrt(15))/180 (10+3*sqrt(15))/45  5/36];

c = [(5-sqrt(15))/10 1/2 (5+sqrt(15))/10];

b   = [5/18 4/9 5/18];
be  = [-5/6 8/3 -5/6];
    
omega = 1.0;

N = length(y0);

I = eye(N);

Jg = @(tn,yn,hn,knsum)[ hn*a(1,1)*(Jf(tn + c(1) * hn, yn + hn * knsum(1,:)))-I,    hn*a(1,2)*(Jf(tn + c(2) * hn, yn + hn * knsum(1,:))),           hn*a(1,3)*(Jf(tn + c(3) * hn, yn + hn * knsum(1,:)));
                        hn*a(2,1)*(Jf(tn + c(1) * hn, yn + hn * knsum(2,:))),           hn*a(2,2)*(Jf(tn + c(2) * hn, yn + hn * knsum(2,:)))-I,    hn*a(2,3)*(Jf(tn + c(3) * hn, yn + hn * knsum(2,:)));
                        hn*a(3,1)*(Jf(tn + c(1) * hn, yn + hn * knsum(3,:))),           hn*a(3,2)*(Jf(tn + c(2) * hn, yn + hn * knsum(3,:))),           hn*a(3,3)*(Jf(tn + c(3) * hn, yn + hn * knsum(3,:)))-I];

                   
while t(n)<dT(2)
    
    % estimation of initial k
    k(1,:) = f(t(n), y_out(n,:)).';
    k(2,:) = f(t(n) + c(2) * h, y_out(n,:) + h * a_exp(2,1) * k(1,:)).';    
    k(3,:) = f(t(n) + c(3) * h, y_out(n,:) + h * (a_exp(3,1) * k(1,:) + a_exp(3,2) * k(2,:))).'; 
    
    err = 1;
    while max(err)>1e-12 % Abbruchbedinung einbauen wenn err nicht konvergiert -> kleineres h
    ksum = [(a(1,1) * k(1,:) + a(1,2) * k(2,:) + a(1,3) * k(3,:));
            (a(2,1) * k(1,:) + a(2,2) * k(2,:) + a(2,3) * k(3,:));
            (a(3,1) * k(1,:) + a(3,2) * k(2,:) + a(3,3) * k(3,:))];
    
    g = [f(t(n) + c(1) * h, y_out(n,:) + h * ksum(1,:))-k(1,:).';
         f(t(n) + c(2) * h, y_out(n,:) + h * ksum(2,:))-k(2,:).';
         f(t(n) + c(3) * h, y_out(n,:) + h * ksum(3,:))-k(3,:).'];
      
    Jgn = Jg(t(n),y_out(n,:),h,ksum);
    
    dk_srt = linsolve(Jgn,-1*g);
    dk = [dk_srt(1:N).';
          dk_srt(N+1:2*N).';
          dk_srt(2*N+1:3*N).'];
    
    k = k + omega * dk;
    
    err = abs(omega*h*[c(1)*norm(dk(1,:)); c(2)*norm(dk(2,:)); c(3)*norm(dk(3,:))]);
    end
    
    y_step   = y_out(n,:) + h * (b(1)*k(1,:) + b(2)*k(2,:) + b(3)*k(3,:));
    ye_step  = y_out(n,:) + h * (be(1)*k(1,:) + be(2)*k(2,:) + be(3)*k(3,:));
     
    EST = norm(y_step - ye_step);
    
    if ((EST/h<=tau) || (h<=h_min))
        n = n+1;
        t(n) = t(n-1) + h;
                
        y_out(n,:) = y_step;
                
        h = max(h_min, min(mu*h, rho*((tau/EST)*h^(p+1) )^(1/p)));
               
        if ((t(n)+h)>dT(2))   h=dT(2) - t(n); end

    else
        rejects = rejects + 1;
        h = max(h_min, 0.5*h);
    end  
end
el_RKnia = toc;
disp(['Steps of RK_imp6_auto:       ' num2str(length(t))])
disp(['Time of RK_imp6_auto:        ' num2str(el_RKnia)])
disp(['No. rejects:                 ' num2str(rejects)])
disp(' ')
t=t.';
end
function [t,y_out] = RK_xorder(f, dT, y0, h_start, h_min, tau)
% RK mit anpassung der Ordnung je nach Fehler aus:
% A Variable Order Runge-Kutta Method for Initial Value Problems with Rapidly Varying Right-Hand Sides 
% von Cash

tic

rejects = 0;
n = 1;
t(n) = dT(1);
h = h_start;
y_out(n,:) = y0(:);


twiddle(n,1:2) = [1.5 1.1];
quit(n,1:2) = [100 100];
ERR(n,1:5) = 0;
E(n,1:5) = 0;
SF = 0.9;

n = n+1;

% RK coefficients
% a = [0,             0,          0,          0,              0,        0;...
%      1/5,           0,          0,          0,              0,        0;...
%      3/40,          9/40,       0,          0,              0,        0;...
%      3/10,         -9/10,       6/5,        0,              0,        0;...
%      -11/54,        5/2,       -70/27,      35/27,          0,        0;...
%      1631/55296,   175/512,    575/13824,  44275/110592,   253/4096, 0];
%  
% c = [0; 1/5; 3/10; 3/5; 1; 7/8];
% 
% b1 = [1,            0,      0,          0,              0,          0];
% b2 = [-3/2,         5/2,    0,          0,              0,          0];
% b3 = [19/54,        0,   -10/27,        55/54,          0,          0];
% b4 = [2825/27648,   0,    18575/48384,  13525/55296,    277/14336,  1/4];
% b5 = [37/378,       0,    250/621,      125/594,        0,          512/1771];


while t(n-1)<dT(2)
    k(1,:) = f(t(n-1), y_out(n-1,:)).';
    k(2,:) = f(t(n-1) + 1/5 * h, y_out(n-1,:) + h * 1/5 * k(1,:)).';
    
    y(1,:) = y_out(n-1,:) + h * k(1,:);
    y(2,:) = y_out(n-1,:) + h * (-3/2 * k(1,:) + 5/2 * k(2,:));
    
    ERR(n,1)   = norm(y(2,:) - y(1,:))^(1/2);
    E(n,1)     = ERR(n,1) / tau^(1/2);
    
    if E(n,1) > twiddle(n-1,1) * quit(n-1,1) && h > h_min
        % abandon the step and adjust h
        esttol = E(n,1) / quit(n-1,1);
        h = max(1/5, SF/esttol) * h;
        rejects = rejects + 1;
    else
       k(3,:) = f(t(n-1) + 3/10 * h, y_out(n-1,:) + h * (3/40 * k(1,:) + 9/40 * k(2,:))).';
       k(4,:) = f(t(n-1) + 3/5  * h, y_out(n-1,:) + h * (3/10 * k(1,:) +-9/10 * k(2,:) + 6/5 *k(3,:))).';

       y(3,:) = y_out(n-1,:) + h * (19/54 * k(1,:) + -10/27 * k(3,:) + 55/54 * k(4,:));
       
       ERR(n,2) = norm(y(3,:) - y(2,:))^(1/3);
       E(n,2)   = ERR(n,2) / tau^(1/3);
       
       if E(n,2) > twiddle(n-1,2) * quit(n-1,2)
           % Try a lower order Solution
           if E(n,1)<1
               % check the error of the second order solution
               if norm(h/10 * (k(2,:)-k(1,:)))<tau
                   % accept second order solution
                   y_out(n,:) = y_out(n-1,:) + h/10 * (k(1,:) + k(2,:));
                   h = h/5;
                   t(n) = t(n-1) + h;
                   twiddle(n,:) = twiddle(n-1,:);
                   quit(n,:) = quit(n-1,:);
                   n = n+1;
                   else
                   % abandon the step
                   h = h/5;
                   rejects = rejects + 1;
               end
           else
               % abandon the step
               esttol = E(n,2) / quit(n-1,2);
               h = max(1/5, SF/esttol) * h;
               rejects = rejects + 1;
           end
       else
           k(5,:) = f(t(n-1) +        h, y_out(n-1,:) + h * (-11/54       * k(1,:) + 5/2        * k(2,:) +-70/27      * k(3,:) + 35/27          * k(4,:))).';
           k(6,:) = f(t(n-1) + 7/8  * h, y_out(n-1,:) + h * (1631/55296   * k(1,:) + 175/512	* k(2,:) + 575/13824  * k(3,:) + 44275/110592   * k(4,:) + 253/4096 * k(5,:))).';
          
           y(4,:) = y_out(n-1,:) + h * (2825/27648    * k(1,:) + 18575/48384    * k(3,:) + 13525/55296	* k(4,:) + 277/14336 * k(5,:) + 1/4         * k(6,:));
           y(5,:) = y_out(n-1,:) + h * (37/378        * k(1,:) + 250/621        * k(3,:) + 125/594        * k(4,:) +                    512/1771    * k(6,:));
           
           ERR(n,4) = norm(y(5,:) - y(4,:))^(1/5);
           E(n,4)   = ERR(n,4) / tau^(1/5);
           
           if E(n,4)>1
               % readjust the twiddle factors
               for i=1:2
               if E(n,i)/quit(n-1,i)<twiddle(n-1,i)
                   twiddle(n,i) = max(1.1, E(n,i)/quit(n-1,i));
               else
                   twiddle(n,i) = twiddle(n-1,i);
               end
               end
               if E(n,2)<1
                   % check the accuracy of the third order solution
                   if norm(h/10 * (k(1,:) - 2*k(3,:) + k(4,:)))<tau
                       % accept the order 3 solution
                       y_out(n,:) = y_out(n-1,:) + h * (1/10 * k(1,:) + 2/5 * k(3,:) + 1/10 * k(4,:));
                       h = 3*h/5;
                       t(n) = t(n-1) + h;
                       twiddle(n,:) = twiddle(n-1,:);
                       quit(n,:) = quit(n-1,:);
                   end
               else
                   % try a lower order solution
                   if E(n,1)<1
                       % check the accuracy of the second order solution
                       if norm(h/10 * (k(2,:)-k(1,:)))<tau
                           % accept second order solution
                           y_out(n,:) = y_out(n-1,:) + h/10 * (k(1,:) + k(2,:));
                           h = h/5;
                           t(n) = t(n-1) + h;
                           twiddle(n,:) = twiddle(n-1,:);
                           quit(n,:) = quit(n-1,:);
                           n = n+1;
                       else
                           % abandon the step
                           h = h/5;
                           rejects = rejects + 1;
                       end
                   else
                       % abandon the current step
                        esttol = E(n,4);
                        h = max(1/5, SF/esttol) * h;
                        rejects = rejects + 1;
                   end
               end
           else
               % accept the order 5 solution
               y_out(n,:) = y(5,:);
               t(n) = t(n-1) + h;
               h = min(5, SF/E(n,4)) * h;
                              
               Q(1) = E(n,1)/E(n,4);
               Q(2) = E(n,2)/E(n,4);
               
               for j=1:2
                   if Q(j)>quit(n-1,j)
                       Q(j) = min(Q(j), 10*quit(n-1,j));
                   else
                       Q(j) = max(Q(j), 2/3 * quit(n-1,j));
                   end
                   quit(n,j) = max(1.0, min(10000,Q(j)));
                   twiddle(n,j) = twiddle(n-1,j);
               end
               n = n+1;
           end
       end
    end
    if ((t(n-1)+h)>dT(2))   h=(dT(2) - t(n-1)); end
end
el_ode = toc;
disp(['Steps of RK_xorder:          ' num2str(length(t))])
disp(['Time of RK_xorder:           ' num2str(el_ode)])
disp(['No. rejects:                 ' num2str(rejects)])
disp(' ')
t=t.';
end
function [t,y] = RK_DP(f, dT, y0, h_start, h_min, tau)
% Eingebettetes Runge Kutta Verfahren von Dormand & Prince
% Fehlerschätzungsformel aus Beispiel 2.28 aus:
% http://www.asc.tuwien.ac.at/~melenk/teach/num_DGL_SS08/ode_teil4.pdf

tic

rejects = 0;
mu = 2;
rho = 0.8;
t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);
est_out(t_count) = 0;

p = 5; % RK of 5th order
a = [0,         0,          0,          0,          0,              0,      0;...
     1/5,       0,          0,          0,          0,              0,      0;...
     3/40,      9/40,       0,          0,          0,              0,      0;...
     44/45,     -56/15,     32/9,       0,          0,              0,      0;...
     19372/6561,-25360/2187,64448/6561, -212/729,   0,              0,      0;...
     9017/3168, -355/33,    46732/5247, 49/176,     -5103/18656,    0,      0;...
     35/384,    0,          500/1113,   125/192,    -2187/6784,     11/84,  0];

b4 = [      35/384,     0,      500/1113,    125/192,       -2187/6784,     11/84,      0];
b5 = [  5179/57600,     0,    7571/16695,    393/640,    -92097/339200,  187/2100,   1/40];

c  = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];

while t(t_count)<dT(2)    
    y4  = RK_step(a,b4,c,f,y(t_count,:),t(t_count),h);
    y5  = RK_step(a,b5,c,f,y(t_count,:),t(t_count),h);
 
    EST = norm(y5 - y4);
    
    if ((EST/h<=tau) || (h<=h_min))
        t_count = t_count + 1;
        t(t_count) = t(t_count-1)+h;
        est_out(t_count) = EST;
        y(t_count,:)=y5;
                
        h = max(h_min, min(mu*h, rho*((tau/EST)*h^(p+1) )^(1/p)));
               
        if ((t(t_count)+h)>dT(2))   h=dT(2) - t(t_count); end

    else
        rejects = rejects + 1;
        h = max(h_min, 0.5*h);
    end
end
el_ode = toc;
disp(['Steps of RK_DP:              ' num2str(length(t))])
disp(['Time of RK_DP:               ' num2str(el_ode)])
disp(['No. rejects:                 ' num2str(rejects)])
disp(' ')
t=t.';
end
function [t,y] = RK4_auto(f, dT, y0, h_start, h_min, tau)
% Klassisches Runge Kutta Verfahren mit
% adaptiver Schrittweitensteuerung mit Extrapolation
% Fehlerschätzung aus Algorithmus 2.23 von:
% http://www.asc.tuwien.ac.at/~melenk/teach/num_DGL_SS08/ode_teil4.pdf
tic

rejects = 0;
p = 4; % RK4 vierter Ordnung
mu = 2;
rho = 0.8;
t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);

% RK coefficients
a = [ 0,  0, 0, 0;...
    1/2,  0, 0, 0;...
      0,1/2, 0, 0;...
      0,  0, 1, 0];
b = [1/6, 1/3, 1/3, 1/6];
c = [0; 1/2; 1/2; 1];
    
while t(t_count)<dT(2)
    h2 = 2*h;
    
    y2      = RK_step(a,b,c,f,y(t_count,:),t(t_count),h2);
    y12     = RK_step(a,b,c,f,y(t_count,:),t(t_count),h);
    y_top   = RK_step(a,b,c,f,y12,t(t_count) + h,h);
        
    EST = norm(y2 - y_top)/(1 - 2^(-p));
 
    if (((EST/h2)<=tau) || h<=h_min)      % %(((EST/h2)<=tau) || h<h_min)  Genauigkeit oder h_min erreicht
        est_out(t_count) = EST;
        t_count = t_count + 1;
        t(t_count) = t(t_count-1)+h;
        y(t_count,:)=y12;
        
        est_out(t_count) = EST;
        t_count = t_count + 1;
        t(t_count) = t(t_count-1)+h;
        y(t_count,:)=y_top;
        
        
        h = 0.5 * max(2*h_min, min(mu*h2, rho*((tau/EST)*h2^(p+1) )^(1/p)));
               
        if ((t(t_count)+2*h)>dT(2))   h=(dT(2) - t(t_count)) / 2; end

    else
        rejects = rejects + 1;
        h = max(h_min, 0.5*h);
    end
end
est_out(t_count) = EST;
el_ode = toc;
disp(['Steps of RK4_auto2:          ' num2str(length(t))])
disp(['Time of RK4_auto2:           ' num2str(el_ode)])
disp(['No. rejects:                 ' num2str(rejects)])
disp(' ')
t=t.';
end
function [t,y] = RK4(f, dT, y0, h_start, h_min, tau)
% Klassisches Runge-Kutta verfahren (m=4)
tic

t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);

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
end
el_ode = toc;
disp(['Steps of RK4:                ' num2str(length(t))])
disp(['Time of RK4:                 ' num2str(el_ode)])
disp(' ')
t=t.';
end
function [t,y] = FWE_auto(f, dT, y0, h_start, h_min, tau)
% Forward Euler method mit automatischer Schrittweitenanpassung nach
% Algorithmus 17.12 auf Seite 432 von Numerik-Algorithmen
tic

t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);

y_sym = sym('y_sym',[length(y0) 1]);
syms time;
J = jacobian(f(time,y_sym),y_sym);
f2 = matlabFunction(diff(f(time,y_sym),time) + J * f(time,y_sym),...
    'Optimize',false,...
    'Vars', [time, y_sym(:).']);

iteration = 0;
while t(t_count)<dT(2)
    y_args = num2cell(y(t_count,:));
    
    y1 = y(t_count,:) + h * f(t(t_count),y(t_count,:))';
    y2 = y(t_count,:) + h * f(t(t_count),y(t_count,:)).' + 0.5 * h^2 * f2(t(t_count),y_args{:}).'; %(t_count,y(t_count,1),y(t_count,2),y(t_count,3)).';
    
    S = (h*tau / norm(y1 - y2))^(1/2) ;%
    
    if S>=1||h<=h_min
        y(t_count+1,:) = y2;
               
        t_count = t_count + 1 ;%
        t(t_count) = t(t_count-1)+h ;%
        
        h = max(min(2,S)*h,h_min) ;%
    else
        iteration = iteration+1;
        h = max(max(1/2,S)*h,h_min) ;%
    end
end
el_ode = toc;
disp(['Steps of FWE_auto:           ' num2str(length(t))])
disp(['Time of FWE_auto:            ' num2str(el_ode)])
disp(['FWE_auto Iterations:         ' num2str(iteration)])
disp(' ')
t=t.';
end
function [t,y] = FWE2(f, dT, y0, h_start, h_min, tau)
%Forward Euler Method zweiter Ordnung
tic

t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);

y_sym = sym('y_sym',[length(y0) 1]);
syms time;
J = jacobian(f(time,y_sym),y_sym);
f2 = matlabFunction(diff(f(time,y_sym),time) + J * f(time,y_sym),...
    'Optimize',false,...
    'Vars', [time, y_sym(:).']);


while t(t_count)<dT(2)
    y_args = num2cell(y(t_count,:));
    y(t_count+1,:) = y(t_count,:) + h * f(t(t_count),y(t_count,:)).' + 0.5 * h^2 * f2(t(t_count),y_args{:}).'; %(t_count,y(t_count,1),y(t_count,2),y(t_count,3)).';
    
    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
end
el_ode = toc;
disp(['Steps of FWE2:               ' num2str(length(t))])
disp(['Time of FWE2:                ' num2str(el_ode)])
disp(' ')
t=t.';
end
function [t,y] = FWE(f, dT, y0, h_start, h_min, tau)
% Forward Euler Method
tic
t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);

while t(t_count)<dT(2)
    y(t_count+1,:) = y(t_count,:) + h * f(t(t_count),y(t_count,:))';

    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
end
el_ode = toc;
disp(['Steps of FWE:                ' num2str(length(t))])
disp(['Time of FWE:                 ' num2str(el_ode)])
disp(' ')
t=t.';
end
function [t,y] = BWE_fp(f, dT, y0, h_start, h_min, tau)
% Backward Euler Method with fixed point iteration
tic
h = h_start;
t_count = 1;
t(t_count) = dT(1);
y(t_count,:) = y0(:);
counter = 0;

while t(t_count)<dT(2)
    
yk1 = y(t_count,:);
dy = 100;

    while dy>1e-7&counter<=5
    counter = counter + 1;
    prev = yk1;
    yk1 = y(t_count,:) + h * f(t(t_count) + h, prev).';
    dy = abs(prev - yk1);    
    end
    counter=0;
    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
    
    y(t_count,:) = yk1;
end
el_ode = toc;
disp(['Steps of BWE_fp:             ' num2str(length(t))])
disp(['Time of BWE_fp:              ' num2str(el_ode)])
disp(' ')
t=t.';
end

% Weitere Solver Funktionen
function fig = plotNumeric(sel, ode_case, out, n)
fig = figure('Position', [-788 -179 789 984]);
hold on

if ode_case==4
    % For Orbiting:
    for i=1:length(sel)
        eval(['plot(out.' char(sel(i)) '.yout(:,1), out.' char(sel(i)) '.yout(:,3), ''DisplayName'', char(sel(i)));']);
    end
else
    for i=1:length(sel)
        for it=1:n
        subplot(n,1,it)
        hold on
        
        eval(['plot(out.' char(sel(i)) '.time, out.' char(sel(i)) '.yout(:,it), ''DisplayName'', char(sel(i)));']);
        legend;
        end
   end
end
end
function out = runSolver(sel, f, dT, y0, h_start, h_min, tau)
out = struct;

for i=1:length(sel)
eval(['[t_out, y_out] = ' char(sel(i)) '(f, dT, y0, h_start, h_min, tau);']);  
eval([char(sel(i)) '= struct(''time'',t_out, ''yout'',y_out);']);
[out(:).(char(sel(i)))] = eval(char(sel(i)));
end
end
function [yn1] = RK_step(a,b,c,f,yn,t,h)
    s = length(c);
    k(1:s,1:length(yn)) = 0;
    for i=1:s
        k(i,:)   = f(t + c(i) * h, yn + h * (a(i,:)*k)).';
    end
    yn1 = yn + h * b * k;    
end 
function Y = ode4(odefun,tspan,y0,varargin)
%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);

Y(:,1) = y0;
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  yi = Y(:,i-1);
  F(:,1) = feval(odefun,ti,yi,varargin{:});
  F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),varargin{:});
  F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),varargin{:});  
  F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),varargin{:});
  Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
end
Y = Y.';
end
function [fct, fct_sol, y0] = getfct(n,lbd)
switch n    
    case 1
    %% dahlquistsche Testgleichung
    lmbd = -3;
    fct = @(t,y)[lmbd*y];
    fct_sol = @(t)[exp(lmbd.*t)];
    y0 = ones(1,1);
    
    case 2
    %% stiff ODE system
    fct = @(t,y)[98 * y(1) + 198 * y(2); -99 * y(1) - 199 * y(2)];
    fct_sol = @(t)[2 .* exp(-t) - exp(-100 .* t); -exp(-t) + exp(-100.*t)];
    y0 = [1 0];
    
    case 3
    %% Halogen ODE
    % parameters for a simplyfied halogen lamp system
    R_c = 0.37;      % electrical resistance @ T_u in ohm
    L = 15e-6;       % inductance in H
    c_w = 7e-3;      % heat capacity of filament in J/K
    T_u = 295;       % ambient temperature in K
    b = 1.5e-12;     % coefficient to compute radiant power in W/K^4
    k = 0.9;         % exponent to compute electrical resistance at hot temperatures
    % parameters for general equation form
    %p = zeros(1,6)
    p(1) = -R_c / (L * T_u * k);
    p(2) = k;
    p(3) = 1 / L;
    p(4) = -b / c_w;
    p(5) = T_u;
    p(6) = R_c / (c_w * T_u * k);
    % define DGLs
    fct = @(t,x)[...
        p(1) * x(2) *p(2) * x(1);
        p(4) * (x(2) - p(5)) *4 + p(6) * x(2) * p(2) * x(1) *2];
    y0 = [0; 100];
    fct_sol = 0;

    case 4
    %% Orbiting
    mu = 0.012277471;
    eta = 1 - mu;
    
    % a = (((y(1)+mu)^2+y(3)^2)^(-1/2))
    % b = (((y(1)-eta)^2+y(3)^2)^(-1/2))
    
    fct = @(t,y)[   y(2);
                    y(1)*(1-eta*(((y(1)+mu)^2+y(3)^2)^(-1/2))^3-mu*(((y(1)-eta)^2+y(3)^2)^(-1/2))^3)+2*y(4)-eta*mu*((((y(1)+mu)^2+y(3)^2)^(-1/2))^3-(((y(1)-eta)^2+y(3)^2)^(-1/2))^3);
                    y(4);
                    y(3)*(1-eta*(((y(1)+mu)^2+y(3)^2)^(-1/2))^3-mu*(((y(1)-eta)^2+y(3)^2)^(-1/2))^3)-2*y(2)];
    y0 = [0.994 0 0 -2.00158510637908252240537862224];
    fct_sol = 0;

   case 5
    %% stiff ODE
    fct = @(t,y)[-1000 * y + 3000 - 2000 * exp(-t)];
    fct_sol = @(t)[3 - 0.998 .* exp(-1000.*t) - 2.002 .* exp(-t)];
    y0 = zeros(1);
    
    case 6
    fct = @(t,y)[lbd * y + (1-lbd)*cos(t) - (1+lbd)*sin(t)];
    fct_sol = @(t)[sin(t) + cos(t)];
    y0 = ones(1);

    case 7
    fct = @(t,y)[-2*y(1)^2;    -3*y(2)^3;     -4*y(3)^2];
    fct_sol = @(t)[(1./(2.*t + 1)); ((2.^(1/2)*(1./(3*t + 1/2)).^(1/2))./2); (1./(4.*t + 1))];
    y0 = ones(3,1);

    case 8
    fct = @(t,y)[-200*t*y^2];
    fct_sol = @(t)[1./(1+100.*t.^2)];
    y0 = ones(1,1);
end
end
function [err] = abserr(f, t, y)
f_val = f(t.').';
err = abs(f_val - y);
end

% ODE analyzation functions
function [eigV] = comp_eigenvalues(f, r, n)
    % compute eigenvalues of non linear DGL system"""

    syms t;
    x = sym('x',[1 size(r,1)]);
    Jf = jacobian(sym(f(t, x)),x);
    %Jf = matlabFunction(Jf,'Vars',x);
    xe = gen_xe(n,r);
    N = size(r,1);

    % solve A for all combinations of xe entrys and compute eigenvalues
    for i=1:size(xe,1)
        A_temp = Jf;
        for j=1:N
            A_temp = subs(A_temp, x(j), xe(i,j));
        end
        eigV(i,:) = eig(double(A_temp));
    end
end
function xe = gen_xe(n,r)
    %generate array with all combinations of x[i]
    for i=1:size(r,1)
        xe(i,:) = linspace(r(i,1), r(i,2), n);
    end
    xe = cartprod(xe);
end 
function X = cartprod(matrixIn)
%CARTPROD Cartesian product of multiple sets.
%
%   X = CARTPROD(A,B,C,...) returns the cartesian product of the sets 
%   A,B,C, etc, where A,B,C, are numerical vectors.  
%
%   Example: A = [-1 -3 -5];   B = [10 11];   C = [0 1];
% 
%   X = cartprod(A,B,C)
%   X =
% 
%     -5    10     0
%     -3    10     0
%     -1    10     0
%     -5    11     0
%     -3    11     0
%     -1    11     0
%     -5    10     1
%     -3    10     1
%     -1    10     1
%     -5    11     1
%     -3    11     1
%     -1    11     1
%
%   This function requires IND2SUBVECT, also available (I hope) on the MathWorks 
%   File Exchange site.

varargin = num2cell(matrixIn.',1);

numSets = size(matrixIn,1);
for i = 1:numSets
    thisSet = sort(varargin{i});
    if ~isequal(prod(size(thisSet)),length(thisSet))
        error('All inputs must be vectors.')
    end
    if ~isnumeric(thisSet)
        error('All inputs must be numeric.')
    end
    if ~isequal(thisSet,unique(thisSet))
        error(['Input set' ' ' num2str(i) ' ' 'contains duplicated elements.'])
    end
    sizeThisSet(i) = length(thisSet);
    varargin{i} = thisSet;
end

X = zeros(prod(sizeThisSet),numSets);
for i = 1:size(X,1)
    
    % Envision imaginary n-d array with dimension "sizeThisSet" ...
    % = length(varargin{1}) x length(varargin{2}) x ...
    
    ixVect = ind2subVect(sizeThisSet,i);
    
    for j = 1:numSets
        X(i,j) = varargin{j}(ixVect(j));
    end
end
end
function X = ind2subVect(siz,ndx)
%IND2SUBVECT Multiple subscripts from linear index.
%   IND2SUBVECT is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   X = IND2SUBVECT(SIZ,IND) returns the matrix X = [I J] containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%
%   For N-D arrays, X = IND2SUBVECT(SIZ,IND) returns matrix X = [I J K ...]
%   containing the equivalent N-D array subscripts equivalent to IND for 
%   an array of size SIZ.
%
%   See also IND2SUB.  (IND2SUBVECT makes a one-line change to IND2SUB so as
%   to return a vector of N indices rather than retuning N individual
%   variables.)%IND2SUBVECT Multiple subscripts from linear index.
%   IND2SUBVECT is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   X = IND2SUBVECT(SIZ,IND) returns the matrix X = [I J] containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%
%   For N-D arrays, X = IND2SUBVECT(SIZ,IND) returns matrix X = [I J K ...]
%   containing the equivalent N-D array subscripts equivalent to IND for 
%   an array of size SIZ.
%
%   See also IND2SUB.  (IND2SUBVECT makes a one-line change to IND2SUB so as
%   to return a vector of N indices rather than returning N individual
%   variables.)
 

% All MathWorks' code from IND2SUB, except as noted:

n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1
  X(i) = floor(ndx/k(i))+1;      % replaced "varargout{i}" with "X(i)"
  ndx = rem(ndx,k(i));
end
end

% Stability Analysis
function h_max = get_h_max(a,b,c,lmbd)
syms yn h x R;
f(x) = lmbd * x;

    s = length(c);
    k = sym('k', [s 1]);
    for i=1:s
        k(i,:)   = f(yn + h * (a(i,:)*k));
    end
    
   R = simplify((yn + h * b * k)/yn);
   h_sol = double(solve(1 == R, h));
   h_sol = h_sol(imag(h_sol)==0);                % save only real rots
   h_max = max(h_sol);
end

%% Not working
% Function to calculate the analytical equation for xk+1 for implicit solve
% function sol = getXk1(f,n)
%     xk = sym('xk_',[n 1]);
%     xk1 = sym('xk1_',[n 1]);
%     syms tk1 h;
%     eqn = xk1 == xk + h*f(tk1,xk1);
%     for i=1:n
%         soli = vpasolve(eqn(i), xk1(i));
%         sol_temp(i) = soli(1);
%     end
%         sol = matlabFunction(sol_temp);
% end

% Backward euler with xk+1 analytical solution
% function [t,y] = BWE(f, dT, y0, h_start, h_min, tau)
% % Backward Euler Method
% xk1 = getXk1(f,length(y0));
% tic
% h = h_start;
% t_count = 1;
% t(t_count) = dT(1);
% y(t_count,:) = y0(:);
% 
% while t(t_count)<dT(2)
%     y(t_count+1,:) = xk1(t(t_count) + h, y(t_count,:), h);
% 
%     t_count = t_count + 1;
%     t(t_count) = t(t_count-1)+h;
% end
% el_ode = toc;
% disp(['Steps of BWE:           ' num2str(length(t))])
% disp(['Time of BWE:            ' num2str(el_ode)])
% disp(' ')
% t=t.';
% end
