clc;
clear all;
close all;

n = 100;                            % number of space gridpoints
% time step
tf = 2;                                          % final time
method = 3;
a = -3;
b = 3;
x = linspace(a,b,n)';
y = x;
h = x(2)-x(1);
r = 0.99;
dt = r*h;
disp(sprintf('Courant number: %0.2f',r))
u = ic(x,y); u0 = u;u_star = u0;
for tn = 1:ceil(tf/dt)
    c = sign(tn*dt-0.5);
    c = 1;
    switch method
        case 1
            name = 'upwind';
            u_star(2:end,2:end) = u(2:end,2:end) - c*(dt/h)*(u(2:end,2:end) - ...
                u(1:end-1,2:end));
            u = u_star;
            u(2:end,2:end) = u_star(2:end,2:end) - c*(dt/h)*(u_star(2:end,2:end) - ...
                u(2:end,1:end-1));            
        case 2
            name = 'Lax-Friedrichs';
            u_star(2:end-1,2:end-1) = 0.5*(u(3:end,2:end-1) + u(1:end-2,2:end-1)) ...
                - c*0.5*(dt/h)*(u(3:end,2:end-1) - u(1:end-2,2:end-1));
            u = u_star;
            u(2:end-1,2:end-1) = 0.5*(u_star(2:end-1,3:end) + u_star(2:end-1,1:end-2)) ...
                - c*0.5*(dt/h)*(u_star(2:end-1,3:end) - u_star(2:end-1,1:end-2));
        case 3
            name = 'Lax-Wendroff';
            u_star(2:end-1,2:end-1) = u(2:end-1,2:end-1)-c*0.5*(dt/h)*...
                (u(3:end,2:end-1)-u(1:end-2,2:end-1))+ ...
             0.5*(c*dt/h)^2*(u(3:end,2:end-1)-2*u(2:end-1,2:end-1)+u(1:end-2,2:end-1)); 
         u = u_star;  
         u(2:end-1,2:end-1) = u_star(2:end-1,2:end-1)-c*0.5*(dt/h)*...
                (u_star(2:end-1,3:end)-u_star(2:end-1,1:end-2))+ ...
             0.5*(c*dt/h)^2*(u_star(2:end-1,3:end)-2*u_star(2:end-1,2:end-1)+...
             u_star(2:end-1,1:end-2));  
    end
    clf
    [X,Y] = meshgrid(x);
    surf(X,Y,u);
    axis([a b a b -.1 1.1])
    title(sprintf('%s , t=%0.2f',name,tn*dt))
    drawnow
end

function Z = ic(x,y)
% initial condition function
[s,~] = size(x);
[X,Y] = meshgrid(x,y);
pt1 = 0;  %Intervalo da função Indicadora
pt2 = 1;
pos1 = Indicador(x,pt1);
pos2 = Indicador(x,pt2);
x(pos1);
x(pos2);
Z = zeros(s,s);
Z(pos1:pos2,pos1:pos2) = sin(pi*X(pos1:pos2,pos1:pos2)).*...
    sin(pi*Y(pos1:pos2,pos1:pos2)); % função indicadora 1
Z(pos1:pos2,pos1:pos2) = ones((pos2-pos1)+1,(pos2-pos1)+1);
surf(X,Y,Z);
end
function n = Indicador (x,pt)
h = x(2) - x(1);
a = x(1);
n = ceil((pt - a)/h)+1;
end
