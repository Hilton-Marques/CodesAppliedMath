clc;
clear all;
close all;

%Métodos numéricos para u'(t) = f(u,t)
dt = 0.5;
ti = 0;
tf = 1;
n =(tf-ti)/dt;
u0 = 2;
t = zeros(n,1);
u = [u0 u0 u0; zeros(n,1) zeros(n,1) zeros(n,1)];
for i =1:(n+1)
    t(i+1,1) = t(i,1) + dt;
    u(i+1,1) = u(i,1) + dt*f(u(i,2),t(i));
    u(i+1,2) = u(i,2) + dt*(0.5*(f(u(i,2),t(i))+f((u(i,2)+dt*f(u(i,2),t(i))),...
        t(i)+dt)));
    k1 = f(u(i,3),t(i));
    k2 = f(u(i,3)+0.5*dt*k1, t(i) + dt*0.5);
    k3 = f(u(i,3) + 0.5*dt*k2, t(i) + dt*0.5);
    k4 = f(u(i,3) + dt*k3, t(i) + dt);    
    u(i+1,3) = u(i,3) + dt*(1/6)*(k1 + 2*k2 + 2*k3 + k4);
    clf
    plot(t(1:i,1),u(1:i,1),'r.-',t(1:i),u(1:i,2),'b:',t(1:i),u(1:i,3),'g.-')
    title(sprintf('time t=%0.2f',t(i)))
    axis([0 4 1 50]);
    drawnow
    pause(0.1)
end
u(1:n+1,1:3)

function v = f(u,t)
v = u*exp(t) - t - 0.25*u;
end

