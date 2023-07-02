clc;
clear all;
close all;

x = [5,4,3,2,1];
fx = x.^1.2;
hold on
axis equal
plot(ones(1,5),-x,'o');
plot(2*ones(1,5),-fx,'o');
hold off

figure
N  = 20;
t = linspace(0,2,N);
x = linspace(-1,1,10);
for i = 1:N
    y = exp((-x.^2)/(4*t(i)))./(4*pi*t(i))^0.5;
    axis equal
    plot(x,y);
    axis([-1.2 1.2 0.2 1])
    title(sprintf('time t=%0.2f',t(i)))
    drawnow
    pause(0.5);
    clf;
end