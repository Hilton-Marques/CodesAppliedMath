clc;
clear all;
close all;

t = 0.1;
xo = [0;0];
yo = [2;0];
uo = yo - xo ;
step1 = t*[uo(2);-uo(1)];
vector2D(uo,step1+uo,[0,0,1])

Gtf = [1 t/2; -t/2 1];
Gtb = [1 -t/2; t/2 1];
Gt = inv(Gtb)*Gtf; %Trapezoidal matrix

Glb = [1 0;t 1];
Glf = [1 t;0 1];
Gl = inv(Glb)*Glf;

F = [(1-(t^2)/2) t;t*(t^2-4)/4 (1-(t^2)/2)];
A = F;  %Leaprofg tipo Verlet (livro do Feynman)


for i = 1:round(2*pi/t)+1
   un = A*uo
   step = un - uo;
   vector2D(uo,step+uo,[0,0,1])
   uo = un;
   pause(0.05)
end 


function hf = vector2D(xi,xf,c)
x = xf - xi;
X = [xi(1,1), xf(1,1)];
Y = [xi(2,1), xf(2,1)];
hl = line( X, Y, 'Color', c,'Linewidth',2);
axis([-10 10    -10 10])
hold on
ht = triangle(xi,xf,[0,0,1]);
grid on
drawnow
hf = [hl,ht];
end
function ht = triangle(xi,xf,c)
h = 0.25;
b = 0.25;
x = xf - xi;
p_x = [-x(2,1); x(1,1)];
p_x = p_x/norm(p_x);
p1 = xf + b/2*p_x;
p2 = xf - b/2*p_x;
p3 = xf + h*(x/norm(x));

X = [p1(1,1), p2(1,1) , p3(1,1)];
Y = [p1(2,1), p2(2,1) , p3(2,1)];
ht = fill(X, Y, c);
end