clc;
clear;
close all;
figure
hold on
%axis([-1.5,1.5,-1.5,1.5]);
t = 100;
dt = 0.001;
N = t/dt;
p0 = complex(1,0);
for i = 1:N
    p1 = p0 + dt*1i*p0;
    line([real(p0),real(p1)],[imag(p0),imag(p1)]);
    p0 = p1;
    drawnow
    %pause(0.001);
end