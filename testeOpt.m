clc;
clear;
close all;
%% Solução Levy

f= @(x)sin(3.*pi.*x(1)).^2 + (((-1)+x(1)).^2) *( 1 + sin(3.*pi.*x(2)).^2)+(((-1)+x(2)).^2)*(1+ ...
  sin(2.*pi.*x(2)).^2);
x0 = [10,10];
[x,fval] = fminunc(f,x0)
f2 = @(x)-20*exp(-0.2*sqrt(0.5*(x(1)^2 + x(2)^2))) - exp(0.5*(cos(2*pi*x(1)) + cos(2*pi*x(2)))) + exp(1) + 20;
[x,fval] = fminunc(f2,x0)


f2 = @(x,y)(sin(3*pi*x)).^2 + ((x - 1).^2)*(1 + (sin(3*pi*y)).^2)+...
    ((y - 1).^2) * (1 + (sin(2*pi*y)).^2);
[X,Y] = meshgrid(linspace(-10,10,10));
Z = f2(X,Y) ;     
figure
hold on
axis equal
surf(X,Y,Z)
view(30,30)

function levi(x)
n = 2;
for i = 1:n; z(i) = 1+(x(i)-1)/4; end
s = sin(pi*z(1))^2;
for i = 1:n-1
    s = s+(z(i)-1)^2*(1+10*(sin(pi*z(i)+1))^2);
end 
y = s+(z(n)-1)^2*(1+(sin(2*pi*z(n)))^2);
end