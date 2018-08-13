%%
A = [-1 0 0 0 0; 0 -2 0 0 0; 0 0 -0.5 0 0; 0 0 0 -1 0; 0 0 0 0 -1000]; B = [0; 0; 0; 0; 1]; C = [2 1 1 1 1]; D = 0;
x0 = [0.125; 0.05; 0.1; 0.2; 0.3];
SYS = ss(A,B,C,D);
time = 3;

eigen_A = eig(A);
ts = -2/min(real(eigen_A))*0.2;

[t,y,x] = FWE(SYS, ts, time, x0);

% plot(t,y);

figure;
subplot(1,2,1)
hold on
plot(t,y)
for i=1:length(x0)
plot(t,x(i,:));
end

[tb,yb,xb] = BWE(SYS, ts, time, x0);

% hold on
% plot(tb,yb,'--');

subplot(1,2,1)
hold on
plot(tb,yb,'--g')
for i=1:length(x0)
plot(tb,xb(i,:),'--r');
end

subplot(1,2,2)
hold on
plot(complex(eigen_A*ts),'o')
circle(-1,0,1)

axis([min(ts*min(real(eigen_A)), -2) max(ts*max(real(eigen_A)),2) min(ts*min(imag(eigen_A)),-2) max(ts*max(imag(eigen_A)),2)])

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

function [t,y,x] = BWE(SYS, h, T, x0)
N = T/h;
t(1)=0;

for i=1:length(x0)
x(i,1) = x0(i);
end

for n=1:N
x(:,n+1) = inv(eye(length(SYS.A)) - h * SYS.A) * (x(:,n) + h * SYS.B * 1);
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