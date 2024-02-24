clc;
clear all;
close all;


S = struct('Upwind',1, ...
    'LaxF',2,...
    'LaxW',3,...
    'Implicit',4);
method = S.Upwind;         % Select numerical method (1) Upwind (2)LF (3)LW
%Grid Assemble
tf = 2;                                          % final time
a = -3;                                         %initial space
b = 3;                                          %final space
dt = 0.05;                                      %time step
t = linspace(0,tf,tf/dt);
r = 1;                                      %CFL number
c = 1;
h = dt/(max(c)*r); 
n = ceil((b-a)/h);
n2D = n^2;
x = linspace(a,b,n)';
y = x;
h = x(2)-x(1);                                % space step
disp(sprintf('Courant number: %0.2f',r))
%Initial condiiton 
u = ic(x,y); u0 = u;u_star = u0;
%MatrixAssemble
I = speye(n);
I2D = speye(n2D);
e=ones(n,1);
R = diag(ones(1,n-1),-1);
R = sparse(diag(ones(1,n-1),-1));
Dxc = (R-R')/(2*h);
Dxx =(1/h^2)*spdiags([e,-2*e,e],-1:1,n,n);
Dxx2Dx = kron(I,Dxx);
Dxx2Dy = kron(Dxx,I);
Dxc2Dx = kron(I,Dxc);
Dxc2Dy = kron(Dxc,I);
R2Dx = kron(I,R);
R2Dy = kron(R,I);
for tn = 1:ceil(tf/dt)
    c = sign(tn*dt-0.5);
    switch method
        case 1
            name = 'Upwind';
            A = (c >=0)*c*(R2Dx-I2D)/h + ...
                (c < 0)*c*(I2D-R2Dx')/h;
            B = (c >=0)*c*(R2Dy-I2D)/h + ...
                (c<0)*c*(I2D-R2Dy')/h;
            M = (I2D+dt*A);
            u_star = M*u;
            M = (I2D+dt*B);
            u = M*u_star;   
        case 2
            name = 'Lax-Friedrichs';
            A = c*Dxc2Dx+(h^2/dt/2)*Dxx2Dx;
            B = c*Dxc2Dy+(h^2/dt/2)*Dxx2Dy;
            M = (I2D+dt*A);
            u_star = M*u;
            M = (I2D+dt*B);
            u = M*u_star;     
        case 3
            name = 'Lax-Wendroff';
            A = c*Dxc2Dx+(c^2*dt/2)*Dxx2Dx;
            B = c*Dxc2Dy+(c^2*dt/2)*Dxx2Dy;
            M = (I2D+dt*A);
            u_star = M*u;
            M = (I2D+dt*B);
            u = M*u_star; 
        case 4
            name = 'Implicit';
            A = c*Dxc2Dx';
            B = c*Dxc2Dy';
            M = (I2D+dt*A);
            u_star = M\u;
            M = (I2D+dt*B);
            u = M\u_star; 
    end
    clf
    [X,Y] = meshgrid(x);
    Z = griddata(reshape(X',[],1),reshape(Y',[],1),u,X,Y);
    Z0 = griddata(reshape(X',[],1),reshape(Y',[],1),u0,X,Y);
    surf(X,Y,Z);
    title(sprintf('%s , t=%0.2f',name,tn*dt),'FontSize',12)
    axis([-3 3 -3 3 -.1 1.1])
    drawnow
    [posX,posY] = pos(X,Y,tn*dt);
    Zexact = ic2(posX,posY); 
    %surf(X,Y,Zexact);
    %title(sprintf('%s , t=%0.2f','Exact',tn*dt),'FontSize',12)
    %Error
    e = (Zexact - Z);
    Err(tn,:) = norm(e,inf);
end
max(Err)
function Z = ic(x,y)
% initial condition function
[X,Y] = meshgrid(x,y);
Z = ((X >=0 & X <=1) & (Y >=0 & Y <=1)).*sin(pi*X).*...
    sin(pi*Y);
Z = Z(:);
end
function [posX,posY] = pos(X,Y,t)
posX = X - (t>=0.5)*(t - 1) + (t<0.5)*t ;
posY = Y - (t>=0.5)*(t - 1) + (t<0.5)*t ;
end
function Z = ic2(X,Y)
% initial condition function
Z = ((X >=0 & X <=1) & (Y >=0 & Y <=1)).*sin(pi*X).*...
    sin(pi*Y);
end
