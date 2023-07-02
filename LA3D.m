clc;
clear all;
close all;
vi = [0;0;0];
vf = [1;2;3];
v1 = vf - vi;
pause(2)
h1 = vector3D(vi,vf,[1 0 0]);
pause(2)
vi = [0;0;0];
vf = [1;1;1];
v2 = vf - vi;
h1 = vector3D(vi,vf,[1 0 0]);
pause(2)
drawPlan(v1,v2)
vi = [0;0;0];
vf = [1;2;2];
v3 = vf - vi;
pause(2)
vector3D(vi,vf,[1 0 0]);
T = proj(v1,v2,v3);
M = [v1,v2];
x = M\(T*v3);
%e = v3 - T*v3;
%vector3D(T*v3,v3,[1 0 0]);
pause(2)
sumV(x(1)*v1,x(2)*v2)

%Transformation(T,v);

function hg = piramide(xi,xf,color)
scaleH = 8; %Decrease the size height
scaleB = 32; %Decrease the size 's base
x = xf - xi;
a = x(1);
b = x(2);
c = x(3);
vx_plano = -b;
vy_plano = a;
if c ~= 0
    v1_plano = [vx_plano; vy_plano; (a*(-vx_plano+xf(1))+...
        b*(-vy_plano + xf(2))+c*xf(3))/c];
else
    v1_plano =[vx_plano; vy_plano; 0];
end
v2_plano = cross(v1_plano,x);
x_un = x/(vecnorm(x)*scaleH);
v1_plano = v1_plano/(vecnorm(v1_plano)*scaleB);
v2_plano = v2_plano/(vecnorm(v2_plano)*scaleB);

p1 = xf + v1_plano;
p2 = xf - v1_plano;
p3 = xf + v2_plano;
p4 = xf - v2_plano;
p5 = xf + x_un;

%Base
X = [p1(1,1), p3(1,1) , p2(1,1),p4(1,1)];
Y = [p1(2,1), p3(2,1) , p2(2,1),p4(2,1)];
Z = [p1(3,1), p3(3,1) , p2(3,1),p4(3,1)];
h1 = fill3(X, Y,Z, color);
%Lado 1
hold on
X = [p1(1,1), p3(1,1) , p5(1,1)];
Y = [p1(2,1), p3(2,1) , p5(2,1)];
Z = [p1(3,1), p3(3,1) , p5(3,1)];
h2 = fill3(X, Y,Z, color);
%Lado 2
hold on
X = [p3(1,1), p2(1,1) , p5(1,1)];
Y = [p3(2,1), p2(2,1) , p5(2,1)];
Z = [p3(3,1), p2(3,1) , p5(3,1)];
h3 = fill3(X, Y,Z, color);
%lado 3
hold on
X = [p2(1,1), p4(1,1) , p5(1,1)];
Y = [p2(2,1), p4(2,1) , p5(2,1)];
Z = [p2(3,1), p4(3,1) , p5(3,1)];
h4 = fill3(X, Y,Z, color);
hold on
X = [p4(1,1), p1(1,1) , p5(1,1)];
Y = [p4(2,1), p1(2,1) , p5(2,1)];
Z = [p4(3,1), p1(3,1) , p5(3,1)];
h5 = fill3(X, Y,Z, color);
hg = [h1,h2,h3,h4,h5];
end

function ht = vector3D(xi,xf,c)
x = xf - xi;
X = [xi(1,1), xf(1,1)];
Y = [xi(2,1), xf(2,1)];
Z = [xi(3,1), xf(3,1)];
hl = line( X, Y,Z, 'Color', c,'Linewidth',2);
hold on
hp = piramide(xi,xf,[0,0,1]);
axis([-0.5 3    -0.5 3    -1  4])
grid on
view(-155,-66)
drawnow
ht = [hl,hp];
end
%v1,v2 are random vector on the plane
function drawPlan (v1,v2)
color = [ 1 0 0];
scale = 10;
u1 = v1/norm(v1);
projv1 = dot(v1/norm(v1),v2)*(v1/norm(v1));
u2 = (v2-projv1)/norm(v2-projv1);
dot(u1,u2)
p1 = u1*scale;
p2 = -u1*scale;
p3 = u2*scale;
p4 = -u2*scale;
X = [p1(1,1), p3(1,1) , p2(1,1),p4(1,1)];
Y = [p1(2,1), p3(2,1) , p2(2,1),p4(2,1)];
Z = [p1(3,1), p3(3,1) , p2(3,1),p4(3,1)];
h1 = fill3(X, Y,Z, color);
set(h1, 'facealpha',0.5);
end

function Transformation (T,v)
[~,c] = size(v);
nFrames = 100; % number of frames
% Método de Interpolação
collect_Frames = zeros(3,nFrames*c);
count = 1;
for i = 1:nFrames
    R = eye(3) + (i/nFrames)*(T-eye(3));
    nv = R*v;
    collect_Frames(1:3,count:count+c-1) = nv;
    count = count+c;
end
count = 1;
for s = 1:nFrames
    vf = collect_Frames(1:3,count:count+c-1);
    h2 = vector3D([0;0;0],vf,[1 0 0]);
    pause(0.005)
    count = count + c;
    delete(h2);
end
u=T*v;
vector3D([0;0;0],u,[1 0 0]);
end
%v1,v2 belong to the subspace, v3 are out
function T = proj(v1,v2,v3)
u1 = v1/norm(v1);
projv1 = dot(v1/norm(v1),v2)*(v1/norm(v1));
u2 = (v2-projv1)/norm(v2-projv1);
A = [u1,u2];
T = A*A';
Transformation(T,v3)
end
function sumV(v1,v2)
xi = v1;
xf = v1 + v2;

vector3D([0;0;0],v1,[1 1 0]);
pause(1)
vector3D(xi,xf,[1 1 0]);

end