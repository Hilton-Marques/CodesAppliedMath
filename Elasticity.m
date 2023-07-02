clc;
clear all;
close all;

%Object
set(gca,'XColor', 'none','YColor','none','ZColor','none');
O = rec3D(1,1,1); %set largura, comprimento e altura

%Transformation
% T = [2 0 0; 0 -0.5 0; 0 0 -0.5];
%T = [-1 1/20 0; 1/20 -1 0; 0 0 0];
%T = [1 0 0; -0.2 1 -1; -0.2 1 1];
O_unique = unique(O','row','stable')';

u = (O_unique(:,7) + O_unique(:,2))*0.5;
v = (O_unique(:,6) + O_unique(:,4))*0.5;
eixo = u - v;
for i = 1:size(O_unique,2)
    p = O_unique(:,i);
    if dot(eixo,p) == 0
        check = true;
    end
end
hold on
show(O);
%plotAxis(eixo);
vec([0;0;0],eixo)
showDiags(O)
drawPlan(eixo);
axis([-1 1 -1 1 -1 1]);
exportgraphics(gca,'halfplanes.png','Resolution',1000);
T = Rot3D(pi,eixo);
%Polar decomposition
[A,S,B] = svd(T);
U = B*S*B';
R = A*B';
TransformCube(T,O,U,R,1);
%Vector3D
% vi = [0;0;0];
% vf = [0.5;0.5;0.5];
% h1 = vector3D(vi,vf,[1 0 0]);


function O = rec3D(L,w,h)
unitCube = [-0.5 -0.5 -0.5; 0.5 -0.5 -0.5; -0.5 0.5 -0.5;0.5 0.5 -0.5; ...
    -0.5 -0.5 0.5;0.5 -0.5 0.5;0.5 0.5 0.5;-0.5 0.5 0.5];
a = -0.5;
b = 0.5;
n = 1;
x = a:(b-a)/n:b;
y = x;
z = x;
c = b*ones(1,n+1);
fac1 = [[-x;c;c] [-c;c;-z] [x;c;-c] ... 
        [c;c;z]];
fac2 = [[c;-y;c] [c;c;-z] [c;y;-c] ... 
        [c;-c;z]];
fac3 = [[-x;-c;c] [-c;-c;-z] [x;-c;-c] ... 
        [c;-c;z]];
fac4 = [[-c;-y;c] [-c;c;-z] [-c;y;-c] ... 
        [-c;-c;z]];
fac5 = [[-x;-c;c] [-c;-y;c] [x;c;c] ... 
        [c;y;c]];
fac6 = [[-x;-c;-c] [-c;-y;-c] [x;c;-c] ... 
        [c;y;-c]];
unitCube = [fac1 fac2 fac3 fac4 fac5 fac6];
v = 1:4*(n+1); %4 arestas
i = 4*(n+1)*ones(1,4*(n+1));
% fac = [v;v+i;v+2*i;v+3*i;v+4*i;v+5*i];
% h2 = patch('Vertices',unitCube','Faces',fac,...
%     'FaceVertexCData',hsv(6),'FaceColor',[0 0 1]);
% view(-38,41);
I = eye(3);
T = I;
T(1:4:end) = [L w h];
O = T*unitCube;
end
function X = TransformCube(T,O,U,R,flag)
nFrames = 100; % number of frames
%Eule angles
eul = rotm2eul(R);
%Initial position
[~,n]=size(O);
v = 1:(n/6); %6 faces
i = (n/6)*ones(1,(n/6));
fac = [v;v+i;v+2*i;v+3*i;v+4*i;v+5*i];
O_unique = unique(O','row','stable')';
h1 = patch('Vertices',O','Faces',fac,...
     'FaceVertexCData',hsv(6),'FaceColor',[0 0 1]);
 for j = 1:size(O_unique,2)
     p = O_unique(:,j);
     text(p(1),p(2),p(3),num2str(j),'FontSize',12,'VerticalAlignment','top');
 end
view(-38,41);
axis([-1 1 -1 1 -1 1]);
set(h1, 'facealpha',0)
pause(1);
% Inbetween frames
if flag == 1
for i = 1:nFrames
    InStre = eye(3) + (i/nFrames)*(U-eye(3)); %Interpolate Stretch
    InAngle = rotm2eul(eye(3)) + (i/nFrames)*(eul-zeros(1,3));
    InRot = eul2rotm(InAngle);
    nv = InRot*InStre*O;
    h2 = patch('Vertices',nv','Faces',fac,...
     'FaceVertexCData',hsv(6),'FaceColor',[0 0 1]);
    nv_unique = unique(nv','row','stable')';
    for j = 1:size(nv_unique,2)
        p = nv_unique(:,j);
        h2(end+1) = text(p(1),p(2),p(3),num2str(j),'FontSize',15,'Color','red','VerticalAlignment','bottom');
    end
    pause(0.03)
    delete(h2);
end
else
    for i = 1:nFrames
    InStre = eye(3) + (i/nFrames)*(U-eye(3)); %Interpolate Stretch
    nv = InStre*O;
    h2 = patch('Vertices',nv','Faces',fac,...
     'FaceVertexCData',hsv(6),'FaceColor',[0 0 1]);
    pause(0.03)
    delete(h2);
    end
    
    for i = 1:nFrames
    InAngle = rotm2eul(eye(3)) + (i/nFrames)*(eul-rotm2eul(eye(3)));
    InRot = eul2rotm(InAngle);
    I = InRot*nv;
    h2 = patch('Vertices',I','Faces',fac,...
     'FaceVertexCData',hsv(6),'FaceColor',[0 0 1]);
    pause(0.01)
    delete(h2);
    end
end
%Final position
X=T*O;
h2 = patch('Vertices',X','Faces',fac,...
    'FaceVertexCData',hsv(6),'FaceColor',[0 0 1]);
X_unique = unique(X','row','stable')';
for j = 1:size(X_unique,2)
    p = X_unique(:,j);
    h2(end+1) = text(p(1),p(2),p(3),num2str(j),'Color','red','VerticalAlignment','bottom');
end
end
function ht = vector3D(xi,xf,c)
x = xf - xi;
X = [xi(1,1), xf(1,1)];
Y = [xi(2,1), xf(2,1)];
Z = [xi(3,1), xf(3,1)];
hl = line( X, Y,Z, 'Color', c,'Linewidth',2);
hold on
hp = piramide(xi,xf,[0,0,1]);
axis([0 1    0 1    0  1])
grid on
view(-155,-66)
drawnow
ht = [hl,hp];
end
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
function Rot3D = Rot3D(angle, axis)
n = axis / vecnorm(axis); 
chute = [n(3,1);-n(1,1);n(2,1)];
plan1 = (eye(3) - n*n')*chute;
plan1 = plan1/ vecnorm(plan1); 
plan2 = cross(plan1,n);
R = [cos(angle) -sin(angle) 0;...
    sin(angle) cos(angle) 0; ...
    0 0 1];
M = [plan1 plan2 n];
Rot3D = M*R*M';
end
function plotAxis(v)
a = v;
b = -v;
line([a(1),b(1)],[a(2),b(2)],[a(3),b(3)]);
end
function show(O)
[~,n]=size(O);
v = 1:(n/6); %6 faces
i = (n/6)*ones(1,(n/6));
fac = [v;v+i;v+2*i;v+3*i;v+4*i;v+5*i];
O_unique = unique(O','row','stable')';
h1 = patch('Vertices',O','Faces',fac,...
     'FaceVertexCData',hsv(6),'FaceColor',[0 0 1]);
%  for j = 1:size(O_unique,2)
%      p = O_unique(:,j);
%      text(p(1),p(2),p(3),num2str(j),'FontSize',12,'VerticalAlignment','top');
%  end
view(-170,21);
axis([-1 1 -1 1 -1 1]);
set(h1, 'facealpha',0)
end
function showDiags(O)
ids = [7,2,5,1];
for i = 1:4
    p_a = O(:,ids(i));
    plotAxis(p_a);
    text(p_a(1),p_a(2),p_a(3),num2str(i),'FontSize',12,'VerticalAlignment','top');
    text(-p_a(1),-p_a(2),-p_a(3),num2str(i),'FontSize',12,'VerticalAlignment','top');
end
end

function drawPlan(n,d,color,A,xo)
if nargin == 1
    A = 20000;
    xo = [0;0;0];
    color = [1,0,0];
    d = 0;
elseif nargin == 2
    xo = [0;0;0];
    color = [1,0,0];
    A = 20000;
elseif nargin == 3
    xo = [0;0;0];
    A = 2000000;
end
[x,y] = findTriedro(n);
xp1 = xo + A*x - n*d ;
yp1 = xo + A*y - n*d ;
xp2 = xo - A*x - n*d ;
yp2 = xo - A*y - n*d ;
h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.3);
end
function [x,y] = findTriedro(z)
z = z/norm(z);
uTemp = [z(3);z(1);-z(2)];
uTemp = uTemp/norm(uTemp);
if (dot(z,uTemp) == 1)
    uTemp = uTemp([2,1,3]);
end
x = cross(uTemp,z);
y = cross(z,x);
end
function vec(p1,p2)
L = norm(p2 - p1);
axis = p2 - p1;
p_trans = p1 + axis;
line([p1(1),p_trans(1)],[p1(2),p_trans(2)],...
                [p1(3),p_trans(3)],'Color', 'black','Linewidth',1);

axis = axis / norm(axis);
axis_rot = cross(axis,[0,0,1]);
angle =  acos(dot(axis,[0,0,1]));
rot_matrix = axang2rotm([axis_rot,-angle]);


hold on
[X,Y,Z]=cylinder([0 (0.5)],20 );
[m,n] = size(X);
p = -0.05*[X(:)';Y(:)';Z(:)'];
p = rot_matrix * p;
X = reshape(p(1,:),m,n) + p_trans(1);
Y = reshape(p(2,:),m,n) + p_trans(2);
Z = reshape(p(3,:),m,n) + p_trans(3);
h = surf(X,Y,Z,'FaceColor','magenta');
end