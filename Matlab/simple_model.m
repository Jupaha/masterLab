% Parameter vom Feder/Dämpfer system
d = 5;
k = 100;
m = 1;

% Input Kraft u
u = 10;

% Differentialgleichung y(1) = x, y(2) = x_dot
% u = 10 für t<0.1
f = @(t, y)[y(2); (-d/m*y(2) - k/m*y(1) + ((t<0.1)*u))];
% Anfangswerte:
y0 = [0 0];

% Stop Bedingung ist im Löser implementiert ab Zeile 60 (Uncomment)

% Runge Kutta Löser
[tout, yout] = RK_DP(f, [0 10], y0, 1e-3, 1e-5, 1e-7, @intReset);

%Plot
subplot(2,1,1)
plot(tout, yout(:,1))
subplot(2,1,2)
plot(tout, yout(:,2))

function y_reset = intReset(y_in)
y_reset = y_in;
if y_in(1)<=0
    y_reset(2) = 0;
end
end
function [t,y] = RK_DP(f, dT, y0, h_start, h_min, tau, resetFunction)
% Eingebettetes Runge Kutta Verfahren von Dormand & Prince
% Fehlerschätzungsformel aus Beispiel 2.28 aus:
% http://www.asc.tuwien.ac.at/~melenk/teach/num_DGL_SS08/ode_teil4.pdf

disp('RK_DP')
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
    %% Stop Bedingung für x<0
    if exist('resetFunction','var')
        y(t_count,:) = resetFunction(y(t_count,:));
    end
    
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
disp(['-Steps:           ' num2str(length(t))])
disp(['-Time:            ' num2str(el_ode)])
disp(['-Rejects:         ' num2str(rejects)])
disp(' ')
t=t.';
end
function [yn1] = RK_step(a,b,c,f,yn,t,h)
    s = length(c);
    k(1:s,1:length(yn)) = 0;
    for i=1:s
        k(i,:)   = f(t + c(i) * h, yn + h * (a(i,:)*k)).';
    end
    yn1 = yn + h * b * k;    
end 
