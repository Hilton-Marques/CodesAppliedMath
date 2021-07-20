clc;
clear all;
close all;

S = struct('Upwind',1, ...
    'LaxF',2,...
    'LaxW',3,...
    'Implicit',4,...
    'TransientePlot', 1, ...
    'NonTransientePlot',0);
method = S.Upwind;         % Select numerical method (1) Upwind (2)LF (3)LW
%Grid Assemble
tf = 2;                                          % final time
a = -3;                                         %initial space
b = 3;                                          %final space
dt = 0.1;                                      %time step
t = linspace(0,tf,tf/dt);
c = cos(pi*t);                                 %velocity function
r = 1.0;                                      %CFL number
h = dt/(max(c)*r);
x = a:h:b;
n = length(x);                               % space step
disp(sprintf('Courant number: %0.2f',r))
%Iniital Condition

u = ic(x');
%Main

F = @(n,h,dt,u,x,flag)main(n,method,h,x,dt,u,tf,flag);
[U,methodName]=F(n,h,dt,u,x,S.TransientePlot);
%Plot

[X,T] = meshgrid(x,t);
uexact = icExa(pos(X,T));
figure(2),surf(X,T,U);
title(sprintf('%s, CFL = %0.2f',methodName,r))
figure(3),surf(X,T,uexact);
title(sprintf('%s','Exact'))
%Error
k = 7;
Err = zeros(k,1);
for i = 1:k
    dt = 0.05/i;
    c = cos(pi*t);                                 %velocity function                                     %CFL number
    h(i) = dt/(max(c)*r);
    x = a:h(i):b;
    n = length(x);
    t = linspace(0,tf,ceil(tf/dt));
    [X,T] = meshgrid(x,t);
    u = ic(x');
    [uaprox,~] = F(n,h(i),dt,u,x,S.NonTransientePlot);
    uExact = icExa(pos(X,T));
    e = (uaprox - uExact);
    Err(i,1) = max(e(:));
    if (i>1)
        p(i) = (log(Err(i)/Err(i-1))/log(h(i)/h(i-1)));
    end
end
g = sprintf('%0.4f',p(end));
fprintf('Order of convergence: %s\r\n',g)
Err

function [U,methodName] = main(n,method,h,x,dt,u,tf,flag)
%Matrix Assemble
I = speye(n);
e=ones(n,1);
R = sparse(diag(ones(1,n-1),-1));
Dxc = (-R'+R)/(2*h);    %Centered difference
Dxx =(1/h^2)*spdiags([e,-2*e,e],-1:1,n,n);    %Second difference
%Iniital Condition
u0 = u;
for tn = 1:ceil(tf/dt)
    c = cos(pi*tn*dt);
    switch method
        case 1
            methodName = 'Upwind';
            A = (c>=0)*c*(R-I)/h + (c<0)*c*(I-R')/h;
        case 2
            methodName = 'Lax-Friedrichs';
            A = c*Dxc+(h^2/dt/2)*Dxx;
        case 3
            methodName = 'Lax-Wendroff';
            A = c*Dxc+(c^2)*(dt/2)*Dxx;
        case 4
            methodName = 'Implicit';
            A = c*Dxc';
    end
    M = (I+dt*A);
    if (strcmp(methodName, 'Implicit'))
        u = M\u;
    else
        u = M*u;
    end
    U(tn,:) = u;
    if (flag > 0)
        clf
        plot(x,u0,'b:',x,u,'r.-')
        axis([-3 3 -.1 1.1])
        title(sprintf('%s , t=%0.2f',methodName,tn*dt))
        drawnow
    end
end
end
function y = ic(x)
% initial condition function
y = (x>=0&x<=1).*sin(pi*x);
%y = (x>=-0.5&x<=0.5)*1;
end
function pos = pos(X,T)
pos = X - sin(pi*T)/pi;
end
function y = icExa(X)
% initial condition function
y = (X>=0&X<=1).*sin(pi*X);
%y = (X>=-0.5&X<=0.5)*1;
end