% syms te x1 x2 x1dot x2dot;
% u = 1;

% A = [0 1;-12 -7]; B = [0; -3.3]; C = [2 1]; D = 0;
% x0 = [0.125; 0.05];
% 
% x1 = 0.275 - 0.55*exp(-3*te) + 0.4 * exp(-4*te);
% x2 =        1.65 *exp(-3*te) - 1.6 * exp(-4*te);

% x1dot = diff(x1);
% x2dot = diff(x2);
% 
% x = [x1; x2];
% xdot = [x1dot; x2dot];
% 
% xdot - A*x + B*u

%%
% A = [-3]; B = [0]; C = [1]; D = 0;
% x0 = [1];
% SYS = ss(A,B,C,D);
% time = 10;
% ts = 1/3;
%%
A = [0 1;-12 -7]; B = [0; 3.3]; C = [2 1]; D = 0;
x0 = [0.125; 0.05];
SYS = ss(A,B,C,D);
time = 3;
ts = 0.05;

[t,y,x] = FWE(SYS, ts, time, x0);
figure;
subplot(1,2,1)
hold on

plot(t,y)
for i=1:length(x0)
plot(t,x(i,:));
end

subplot(1,2,2)
hold on
plot(complex(eig(A)*ts),'o')
circle(-1,0,1)
axis([-3 2 -2 2])

%%
te = [0:0.01:time];
x1 = 0.275 - 0.55*exp(-3*te) + 0.4 * exp(-4*te);
x2 =        1.65 *exp(-3*te) - 1.6 * exp(-4*te);
xana = [x1; x2];
yana = SYS.C*xana;

subplot(1,2,1)
hold on
plot(te,yana,'--r','LineWidth',2)
plot(te,x1,'--r','LineWidth',2)
plot(te,x2,'--r','LineWidth',2)

%%
function [t,y,x] = FWE(SYS, h, T, x0)
N = T/h;
t(1)=0;

for i=1:length(x0)
x(i,1) = x0(i);
end

for n=1:N
f = SYS.A * x(:,n) + SYS.B * 1;
x(:,n+1) = x(:,n) + h*f;
t(n+1) = t(n)+h;
end

y = SYS.C*x;
end

function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);
end