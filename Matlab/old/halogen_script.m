L = 5;
Tu = 30;
k = 3;
Rc = 10;
b = 5;
Cw = 10;

T = 10;
h = 0.1;
% u as step after 1 sec:
u = [zeros(1,1/h) ones(1,((T/h-1/h)+1))];

%%
Tw(1) = Tu;
i(1) = 0;
x(1) = 0;

for n=1:T/h
f1 = 1/Cw * (u(n)^2/(Rc * (Tw(n)/Tu)^k) - b * (Tw(n) - Tu)^4);

Tw(n+1) = Tw(n) + h * f1;

f2 = (u(n) - (i(n) * Rc * (Tw(n)/Tu)^k))/ L;

i(n+1) = i(n) + h * f2;

x(n+1) = x(n) + h;
end
figure(1)
plot(x,i)
figure(2)
plot(x,Tw)

%%
set_param('halogen','Solver','ode1', 'FixedStep', num2str(h), 'StopTime',num2str(T))
simOut = sim('halogen','SimulationMode','normal',...
            'SaveTime','on','TimeSaveName','tout', ...
            'SaveState','on','StateSaveName','xout',...
            'SaveOutput','on','OutputSaveName','yout',...
            'SaveFormat','Array')
figure(1)
hold on
plot(simOut.tout, simOut.xout(:,1), '--r')

figure(2)
hold on
plot(simOut.tout, simOut.xout(:,2), '--r')
