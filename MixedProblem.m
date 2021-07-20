clear 
clc
close all
%CFL1 for wave problem
%dt/(dx) < 1 => dt < dx
dx = 0.1; 
%CFL2 for heat problem
CFL2 = 0.5; %dt/(dx^2)
dt = CFL2*(dx^2);
%check CFL for wave problem 
CFL1 = dt/dx;
disp(sprintf('Wave CFL: %0.2f',CFL1));
%%Grid
x = -2:dx:2;
t = 0:dt:1;
n = length(x);
%%Matrix Assemble
I = speye(n);
e = ones(n,1);
%Centered difference
Dxc= (1/(2*dx))*spdiags([e,0*e,-e],-1:1,n,n);
%Second difference
Dxx = (1/dx^2)*spdiags([e,-2*e,e],-1:1,n,n);
uo = max(sin(pi*x),0);
u(1,:) = uo; ui = uo';
M = ((Dxx + Dxc)*dt + I);
%Boundary condition
M(1,1) = 1; M(1,2) = 0;
M(n,n) = 1; M(n,n-1) = 0;
for i = 1:length(t) - 1
    clf
    ui = M*ui;
    plot(x,uo,'b:',x,ui,'r.-')
    title(sprintf('time t=%0.2f',i*dt))
    axis([-2 2 0 1])
    drawnow
    u(i+1,:) = ui;
end
figure(1),surf(x,t, u);
title(sprintf('Numerical Solution'));