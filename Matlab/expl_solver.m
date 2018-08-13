%% TODO:
% 
% LESEN: file:///Users/Julius/Downloads/Cash-90.pdf
% file:///Users/Julius/Downloads/Cash-90.pdf
% 
% Daruch RK method mit variable order schreiben
% 
% Anschlie�end mit implizitem RK verfahren auseinander setzen... siehe: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Implicit_methods




%% Numerische L�sung durch verschiedene Solver
% http://www.asc.tuwien.ac.at/~melenk/teach/num_DGL_SS08/ode_teil4.pdf
% Explizites l�sen soweit abgeschlossen...
% Als n�chstes vielleicht analyse betrachten... Dann Wahl von explizit /
% implizit
% RK Ordnung an Fehlerabsch�tzung koppeln -> bei kleinem Fehler kleinere
% Ordnung

% Wichtige Autoren: Dormand and Prince, Verner, Shampine, Cash

close all
clear all
disp(' ')
disp('%%%')

%% stiff ODE
% fct = @(t,y)[-1000 * y + 3000 - 2000 * exp(-t)];
% fct_sol = @(t)[3 - 0.998 .* exp(-1000.*t) - 2.002 .* exp(-t)];
% y0 = zeros(1);

%% stiff ODE system
% fct = @(t,y)[98 * y(1) + 198 * y(2); -99 * y(1) - 199 * y(2)];
% fct_sol = @(t)[2 .* exp(-t) - exp(-100 .* t); -exp(-t) + exp(-100.*t)];
% y0 = [1 0];

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

%% Other ODES

% lmbd = -5;
% fct = @(t,y)[lmbd * y + (1-lmbd)*cos(t) - (1+lmbd)*sin(t)];
% fct_sol = @(t)[sin(t) + cos(t)];
% y0 = ones(1);

% fct = @(t,y)[3*y(1)+y(2);y(2)-y(1)+y(2).^4+y(3).^4;y(2)+y(2).^4+y(3).^4+3];
% y0 = ones(3,1);

% fct = @(t,y)[-2*y(1)^2;    -3*y(2)^3;     -4*y(3)^2];
% fct_sol = @(t)[(1./(2.*t + 1)); ((2.^(1/2)*(1./(3*t + 1/2)).^(1/2))./2); (1./(4.*t + 1))];
% y0 = ones(3,1);

% fct = @(t,y)[-200*t*y^2];
% fct_sol = @(t)[1./(1+100.*t.^2)];
% y0 = ones(1,1);

% fct = @(t,y)[(t<=1/3)*sin(t) + (t>1/3)*sin(1/3-t)];
% y0 = zeros(1,1);
 
% lmbd = -1;
% fct = @(t,y)[lmbd*y];
% fct_sol = @(t)[exp(lmbd.*t)];
% y0 = ones(1,1);

% % parameters for a simplyfied halogen lamp system
%     R_c = 0.37;      % electrical resistance @ T_u in ohm
%     L = 15e-6;       % inductance in H
%     c_w = 7e-3;      % heat capacity of filament in J/K
%     T_u = 295;       % ambient temperature in K
%     b = 1.5e-12;     % coefficient to compute radiant power in W/K^4
%     k = 0.9;         % exponent to compute electrical resistance at hot temperatures
% % parameters for general equation form
%     %p = zeros(1,6)
%     p(1) = -R_c / (L * T_u * k);
%     p(2) = k;
%     p(3) = 1 / L;
%     p(4) = -b / c_w;
%     p(5) = T_u;
%     p(6) = R_c / (c_w * T_u * k);
% % define DGLs
%     fct = @(t,x)[...
%         p(1) * x(2) *p(2) * x(1);
%         p(4) * (x(2) - p(5)) *4 + p(6) * x(2) * p(2) * x(1) *2];
%     y0 = [0; 100];

dT = [0 20];
tol = 1e-7;
h_start = 0.005;
h_min = tol;


%% ODE45
opts = odeset('AbsTol', tol, 'relTol', tol);
tic
[t_ode,y_ode] = ode45(fct, dT, y0, opts);
el_ode = toc;
disp(['Steps of ODE45:      ' num2str(length(t_ode))])
disp(['Time of ODE45:       ' num2str(el_ode)])
disp(' ')


%% Implicit Runge Kutta newton iteration


syms t_sym;
y_sym = sym('y_sym', [1,length(y0)]);
Jf_sym = jacobian(fct(t_sym,y_sym));
Jf = matlabFunction(Jf_sym,'Vars',{t_sym,y_sym});

tic
[t_RKni,y_RKni] = RK_newton(fct, dT, y0, h_start, h_min, tol,Jf);
el_RKni = toc;
disp(['Time of RK impl nwt:' num2str(el_RKni)])
disp(' ')

%% Runge Kutta variable Step
% tic
% [t_RK4a2,y_RK4a2,errest_RK4a2] = RK4_auto2(fct, dT, y0, h_start, h_min, tol);
% el_RK4a2 = toc;
% disp(['Time of RK4a2:       ' num2str(el_RK4a2)])
% disp(' ')

%% Runge Kutta auto order
% tic
% [t_RKxa,y_RKxa,errest_RKxa] = RKx_auto(fct, dT, y0, h_start, h_min , tol);
% el_RKxa = toc;
% disp(['Time of RKx auto:    ' num2str(el_RKxa)])
% disp(' ')
%% Runge Kutta auto cash
% tic
% [t_RKac,y_RKac] = RKx_auto(fct, dT, y0, h_start, h_min, tol);
% el_RKac = toc;
% disp(['Time of RK cash:     ' num2str(el_RKac)])
% disp(' ')
%% Runge Kutta Dormand & Prince
% tic
% [t_RKDP,y_RKDP,errest_RKDP] = RK_DP(fct, dT, y0, h_start, h_min , tol);
% el_RKDP = toc;
% disp(['Time of RKDP:        ' num2str(el_RKDP)])
% disp(' ')

%% Backward Euler
% tic
% Xk1 = getXk1(fct,length(y0));
% [t_BWE,y_BWE] = BWE(Xk1, dT, y0, h_start);
% el_BWE = toc;
% disp(['Time of BWE:       ' num2str(el_BWE)])
% disp(' ')

%% Backward Euler
% tic
% [t_BWEfp,y_BWEfp] = BWE_fp(fct, dT, y0, h_start);
% el_BWEfp = toc;
% disp(['Time of BWE_fp:       ' num2str(el_BWEfp)])
% disp(' ')

%% Forward Euler
% tic
% [t_fwe,y_fwe]      = FWE(fct,dT,y0, h_start);
% el_FWE = toc;
% disp(['Time of FWE:       ' num2str(el_FWE)])
% disp(' ')

%% Other Solvers
%[t_fwe2,y_fwe2]    = FWE2(fct, dT, y0, h_start);
%[t_fwea,y_fwea]    = FWE_auto(fct, dT, y0, h_start, 1e-4);
%[t_RK4,y_RK4]      = RK4(fct, dT, y0, h_start);
%[t_RK4a,y_RK4a]    = RK4_auto(fct, dT, y0, h_start, 1e-12);

%% Plot numerical solutions
plotfig = figure('Position', [-788 -179 789 984]);

% For Orbiting:
plot(y_ode(:,1),y_ode(:,3),'DisplayName','ODE45');
hold on
%plot(y_RKDP(:,1),y_RKDP(:,3),'DisplayName','RKDP');
plot(y_RKni(:,1),y_RKni(:,3),'DisplayName','RKNI');
legend

% for i=1:length(y0)
% subplot(length(y0),1,i)
% hold on
% 
% plot(t_ode,y_ode(:,i),'DisplayName','ODE45');
% % plot(t_RK4a2,y_RK4a2(:,i), 'DisplayName','RK4 auto2');
% %plot(t_RKDP,y_RKDP(:,i), 'DisplayName','RKDP');
% plot(t_RKni,y_RKni(:,i),'DisplayName','RKNI');
% %plot(t_RKac,y_RKac(:,i), 'DisplayName','RKac');
% % plot(t_RKxa,y_RKxa(:,i), 'DisplayName','RKxa');
% %plot(t_fwe,y_fwe(:,i), 'DisplayName','FWE');
% % plot(t_BWEfp,y_BWEfp(:,i), 'DisplayName','BWEfp');
% %plot(t_BWE,y_BWE(:,i), 'DisplayName','BWE');
% %plot(t_fwe2,y_fwe2(:,i), 'DisplayName','FWE2');
% %plot(t_fwea,y_fwea(:,i), 'DisplayName','FWE Auto');
% %plot(t_RK4,y_RK4(:,i), 'DisplayName','RK4');
% %plot(t_RK4a,y_RK4a(:,i), 'DisplayName','RK4 auto');
% 
% legend('Location', 'northeast')
% end

%% Plot Exact Solution
if exist('fct_sol')
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
%     abserr_RK4a2 = abserr(fct_sol, t_RK4a2, y_RK4a2);
%     abserr_RKDP = abserr(fct_sol, t_RKDP, y_RKDP);
% 
%     plot(t_ode, abserr_ode, 'DisplayName', 'Absoluter Fehler ODE')
%     plot(t_RK4a2, abserr_RK4a2, 'DisplayName', 'Absoluter Fehler RK4')
%     plot(t_RKDP, abserr_RKDP, 'DisplayName', 'Absoluter Fehler RKDP')
% end
% 
% legend

%% Error Estimation
% figure(errfig)
% hold on
% plot(t_RK4a2,errest_RK4a2, '--*', 'DisplayName', 'Fehlersch�tzung RK4')
% plot(t_RKDP,errest_RKDP, '--*', 'DisplayName', 'Fehlersch�tzung RKDP')
% 

% %% Step size
% stepfig = figure('Position', [670 0 770 800]);
% 
% t_out = t_RKDP;
% 
% semilogy(t_out(1:end-1), diff(t_out), '-x', 'DisplayName', ['steps: ' num2str(length((t_out)))])
% %axis([0 1 1 1e-8])
% legend
%%
function [yn1] = RK_step(a,b,c,f,yn,t,h)
    s = length(c);
    k(1:s,1:length(yn)) = 0;
    for i=1:s
        k(i,:)   = f(t + c(i) * h, yn + h * (a(i,:)*k)).';
    end
    yn1 = yn + h * b * k;    
end
function [t,y_out] = RK_newton(f, dT, y0, h_start, h_min, eps,Jf)
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
    while max(err)>1e-50
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
disp(['Steps of RK_newton:       ' num2str(length(t))])
%disp(['No. rejects:         ' num2str(rejects)])
t=t.';
end
function [t,y_out] = RKx_auto(f, dT, y0, h_start, h_min, eps)
% RK mit anpassung der Ordnung je nach Fehler aus:
% A Variable Order Runge-Kutta Method for Initial Value Problems with Rapidly Varying Right-Hand Sides 
% von Cash
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
    E(n,1)     = ERR(n,1) / eps^(1/2);
    
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
       E(n,2)   = ERR(n,2) / eps^(1/3);
       
       if E(n,2) > twiddle(n-1,2) * quit(n-1,2)
           % Try a lower order Solution
           if E(n,1)<1
               % check the error of the second order solution
               if norm(h/10 * (k(2,:)-k(1,:)))<eps
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
           E(n,4)   = ERR(n,4) / eps^(1/5);
           
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
                   if norm(h/10 * (k(1,:) - 2*k(3,:) + k(4,:)))<eps
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
                       if norm(h/10 * (k(2,:)-k(1,:)))<eps
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
disp(['Steps of RKac:       ' num2str(length(t))])
disp(['No. rejects:         ' num2str(rejects)])
t=t.';
end
function [t,y,est_out] = RK_DP(f, dT, y0, h_start, h_min , tau)
% Eingebettetes Runge Kutta Verfahren von Dormand & Prince
% Fehlersch�tzungsformel aus Beispiel 2.28 aus:
% http://www.asc.tuwien.ac.at/~melenk/teach/num_DGL_SS08/ode_teil4.pdf
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
disp(['Steps of RKDP:       ' num2str(length(t))])
disp(['No. rejects:         ' num2str(rejects)])
t=t.';
end
function [t,y,est_out] = RK4_auto2(f, dT, y0, h_start, h_min , tau)
% Klassisches Runge Kutta Verfahren mit
% adaptiver Schrittweitensteuerung mit Extrapolation
% Fehlersch�tzung aus Algorithmus 2.23 von:
% http://www.asc.tuwien.ac.at/~melenk/teach/num_DGL_SS08/ode_teil4.pdf
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
disp(['Steps of RK4 auto2:  ' num2str(length(t))])
disp(['No. rejects:         ' num2str(rejects)])
t=t.';
end
function [t,y,est_out] = RK4_auto(f, dT, y0, h_start, h_min, tau)
% Klassisches Runge Kutta Verfahren mit Schrittweitensteuerung nach
% Fehlersch�tzungsformel 17.14 aus numerik-Algorithmen Seite 430-432

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
    y2 = RK_step(a,b,c,f,y(t_count,:),t(t_count),h2);
    
    for i=1:2   
    y(t_count+1,:) = RK_step(a,b,c,f,y(t_count,:),t(t_count),h);
    
    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
    end
    
    err = 1/15 * (y(t_count,:) - y2);
    if abs(err)<=0.9*tau
        h = 1.5*h;%min(1.5*h,0.0025)
        if ((t(t_count)+2*h)>dT(2))   h=(dT(2) - t(t_count)) / 2; end
    elseif abs(err) >= tau
        h = 0.8*h;%max(0.8*h,0.0025)
        if ((t(t_count)+2*h)>dT(2))   h=(dT(2) - t(t_count)) / 2; end
    end
end
disp(['Steps of RK4 auto:    ' num2str(length(t))])
est_out = zeros(length(t));
t=t.';
end
function [t,y] = RK4(f, dT, y0, h_start)
% Klassisches Runge-Kutta verfahren (m=4)
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
disp(['Steps of RK4:       ' num2str(length(t))])
t=t.';
end
function [t,y] = FWE_auto(f, dT, y0, h_start, eps)
% Forward Euler method mit automatischer Schrittweitenanpassung nach
% Algorithmus 17.12 auf Seite 432 von Numerik-Algorithmen

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
    
    S = (h*eps / norm(y1 - y2))^(1/2) ;%
    
    if S>=1
        y(t_count+1,:) = y2;
               
        t_count = t_count + 1 ;%
        t(t_count) = t(t_count-1)+h ;%
        
        h = min(2,S)*h ;%
    else
        iteration = iteration+1;
        h = max(1/2,S)*h ;%
    end
end
disp(['Steps of FWE Auto:   ' num2str(length(t))])
disp(['FWE auto Iterations: ' num2str(iteration)])
t=t.';
end
function [t,y] = FWE2(f, dT, y0, h_start)
%Forward Euler Method zweiter Ordnung
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
disp(['Steps of FWE2:       ' num2str(length(t))])
t=t.';
end
function [t,y] = FWE(f, dT, y0, h_start)
% Forward Euler Method
t_count = 1;
t(t_count) = dT(1);
h = h_start;
y(t_count,:) = y0(:);

while t(t_count)<dT(2)
    y(t_count+1,:) = y(t_count,:) + h * f(t(t_count),y(t_count,:))';

    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
end
disp(['Steps of FWE:        ' num2str(length(t))])
t=t.';
end
function [t,y] = BWE(xk1, dT, y0, h)
% Backward Euler Method
t_count = 1;
t(t_count) = dT(1);
y(t_count,:) = y0(:);

while t(t_count)<dT(2)
    y(t_count+1,:) = xk1(t(t_count) + h, y(t_count,:), h);

    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
end
disp(['Steps of BWE:        ' num2str(length(t))])
t=t.';
end
function [t,y] = BWE_fp(f, dT, y0, h)
% Backward Euler Method
t_count = 1;
t(t_count) = dT(1);
y(t_count,:) = y0(:);
%counter = 0;



while t(t_count)<dT(2)
    
yk1 = y(t_count,:);
dy = 100;

    while dy>1e-7
%     counter = counter + 1;
    prev = yk1;
    yk1 = y(t_count,:) + h * f(t(t_count) + h, prev).';
    dy = abs(prev - yk1);    
    end
    t_count = t_count + 1;
    t(t_count) = t(t_count-1)+h;
    
    y(t_count,:) = yk1;
end
disp(['Steps of BWE_fp:        ' num2str(length(t))])
t=t.';
end
function sol = getXk1(f,n)
    xk = sym('xk_',[n 1]);
    xk1 = sym('xk1_',[n 1]);
    syms tk1 h;
    eqn = xk1 == xk + h*f(tk1,xk1);
    for i=1:n
        soli = vpasolve(eqn(i), xk1(i));
        sol_temp(i) = soli(1);
    end
        sol = matlabFunction(sol_temp);
end
function [err] = abserr(f, t, y)
f_val = f(t.').';
err = abs(f_val - y);
end