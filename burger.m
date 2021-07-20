clc;
clear all;
close all;

S = struct('Upwind',1, ...
    'LaxF',2,...
    'LaxW',3,...
    'Conservative',1,...
    'NonConservative',0,...
    'TransientePlot', 1, ...
    'NonTransientePlot',0);
method = S.Upwind;            % Select numerical method (1) Upwind (2)LF (3)LW
mode = S.Conservative;        % Select mode
%Grid Assemble
tf = 2;                                          % final time
a = -3;                                         %initial space
b = 3;                                          %final space
dt = 0.025;                                      %time step
t = linspace(0,tf,tf/dt);
r = 1;                                      %CFL number
c = 1;
h = dt/(max(c)*r);
n = ceil((b-a)/h);
x = linspace(a,b,n)';
h = x(2)-x(1);                                % space step
disp(sprintf('Courant number: %0.2f',r))
%Initial Condition

u = ic(x);
%Main Function

F = @(h,dt,u,x,flag)main(method,u,dt,h,tf,x,flag,mode);
[U,methodName,modeName]=F(h,dt,u,x,S.TransientePlot);
%Plot

[X,T] = meshgrid(x,t);
uexact = icExa(pos(X,T));
figure(2),surf(X,T,U);
title(sprintf('%s-%s',methodName,modeName))
figure(3),surf(X,T,uexact);
title(sprintf('%s','Exact'))
%Error
k = 3;
Err = zeros(k,1);
for i = 1:k
    dt = 0.05/i;
    h(i) = dt/r;
    x = linspace(a,b,ceil((b-a)/h(i)))';
    t = linspace(0,tf,ceil(tf/dt));
    [X,T] = meshgrid(x,t);
    u = ic(x);
    [uaprox,~,~] = F(h(i),dt,u,x,S.NonTransientePlot);
    uExact = icExa(pos(X,T));
    e = (uaprox - uExact);
    Err(i,1) = norm(e(:));
    if (i>1)
        p(i) = abs(log(Err(i)/Err(i-1))/log(h(i)/h(i-1)));
    end
end

g = sprintf('%0.4f',p(end));
fprintf('Order of convergence: %s\r\n',g)
Err
%Functions
function [U,methodName,modeName] = main(method,u,dt,h,tf,x,flag,mode)
u0 = u;
for tn = 1:ceil(tf/dt)
    switch mode
        case  1
            modeName = 'Conservative';
            switch method
                case 1
                    methodName = 'Upwind';
                    u(2:end) = u(2:end) - (dt/h)*(f(u(2:end)) ...
                        - f(u(1:end-1)));
                case 2
                    methodName = 'Lax-Friedrichs';
                    u(2:end-1) = 0.5*(u(1:end-2) + u(3:end))...
                        - 0.5*(dt/h)*(f(u(3:end)) - f(u(1:end-2)));

                case 3
                    methodName = 'Lax-Wendroff';
                    u(2:end-1) = u(2:end-1) ...
                        - 0.5*dt/h * (f(u(3:end)) - f(u(1:end-2))) ...
                        + 0.5*(dt/h)^2 * ...
                        ( j(0.5*(u(3:end) + u(2:end-1))) .* (f(u(3:end)) - f(u(2:end-1))) - ...
                        j(0.5*(u(2:end-1) + u(1:end-2))) .* (f(u(2:end-1)) - f(u(1:end-2))) );
            end
        case 0
            modeName = 'NonConservative';
            
            switch method              
                
                case 1
                    methodName = 'Upwind';
                    u(2:end-1) = u(2:end-1) - (dt/h)*(u(2:end-1)>=0).*u(2:end-1).*...
                        (u(2:end-1)-u(1:end-2)) - ...
                        (dt/h)*(u(2:end-1)<0).*u(2:end-1).*(u(3:end)-u(2:end-1)) ;
                case 2
                    methodName = 'Lax-Friedrichs';
                    u(2:end-1) = 0.5*(u(3:end)+u(1:end-2)) - (dt/h)*0.5*u(2:end-1).*...
                        (u(3:end)-u(1:end-2));
                case 3
                    methodName = 'Lax-Wendroff';
                    u(2:end-1) = u(2:end-1) - (dt/h)*0.5*u(2:end-1).*...
                        (u(3:end)-u(1:end-2))+ ...
             (dt/h)^2*0.5*(u(2:end-1).^2).*(u(3:end)-2*u(2:end-1)+u(1:end-2));
            end
            
    end
    U(tn,:) = u;
    if (flag>0)
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
function flux = f (u)
flux = 0.5 * u.^2;
end
function jacobian = j(u)
jacobian = u;
end
function pos = pos(X,T)
pos = X - icExa(X).*T;
end
function y = icExa(X)
% initial condition function
y = (X>=0&X<=1).*sin(pi*X);
%y = (X>=-0.5&X<=0.5)*1;
end