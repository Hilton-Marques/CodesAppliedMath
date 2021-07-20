clear
close all
clc
% Choice draw or input
flag = 0;
if flag
data = xlsread('Input','Planilha1');
nPoints = data(1,2);
nEl = data(1,3);
nFaces = nEl;
nNodes = nEl;
Faces(nFaces) = Face();
Elements(nEl) = Element();
Nodes(nEl) = Node();
points = data(2:nPoints+1,1:3);
startEl = nPoints+1;
for i = 1:nFaces
    Faces(i).q = data(startEl+i,4:5);
end
for i = 1:nEl
    Elements(i) = Element(data(startEl+i,1:3),points,Faces(i),i);
end
intNode(1) = Node;
else
nFaces = 6;
nEl = 12;
nNodes = 12;
nIntNodes = 7;
points = [1,0,0;...
    1,0,1;...
    1,1,1;...
    1,1,0;...
    0,0,0;...
    0,0,1;...
    0,1,1;...
    0,1,0];
Faces(nFaces) = Face();
Elements(nEl) = Element();
Nodes(nNodes) = Node();
Faces(1).q = [1;1];
Faces(2).q = [1;2];
Faces(3).q = [1;3];
Faces(4).q = [1;4];
Faces(5).q = [1;5];
Faces(6).q = [1;6];
Elements(1) = Element([1,4,3],points,Faces(1),1);
Elements(2) = Element([1,3,2],points,Faces(1),2);
Elements(3) = Element([4,8,7],points,Faces(2),3);
Elements(4) = Element([4,7,3],points,Faces(2),4);
Elements(5) = Element([8,5,6],points,Faces(3),5);
Elements(6) = Element([8,6,7],points,Faces(3),6);
Elements(7) = Element([5,1,2],points,Faces(4),7);
Elements(8) = Element([5,2,6],points,Faces(4),8);
Elements(9) = Element([4,5,8],points,Faces(5),9);
Elements(10) = Element([4,1,5],points,Faces(5),10);
Elements(11) = Element([6,2,3],points,Faces(6),11);
Elements(12) = Element([6,3,7],points,Faces(6),12);
intNode(nIntNodes) = Node;
intNode(1).pos = [0.5,0.5,1];
intNode(2).pos = [0.5,0.5,0];
intNode(3).pos = [0.5,1,0.5];
intNode(4).pos = [0,0.5,0.5];
intNode(5).pos = [1,0.333,0.6666];
intNode(6).pos = [0.3,0.3,0.7];
intNode(7).pos = [0.5,0.5,0.5];
end

for i = 1:nNodes
    Nodes(i) = Node(Elements(i));
end

obj = Solver(Nodes,Elements,intNode);
tic 
obj.Calculate()
toc
keyboard;
%% Draw Mesh
hold on
for i = 1:nNodes
    p = Nodes(i).pos;
    plot3(p(1),p(2),p(3),'o','Color','Red');
    %text(p(1),p(2),p(3),num2str(i),'HorizontalAlignment','left');
end

for i = 1:nEl
    element = Elements(i);
    p1 = element.points(1,:);
    p2 = element.points(2,:);
    p3 = element.points(3,:);
    plotTriangle(i,p1,p2,p3);
    xc = (p1 + p2 + p3)/3;
    vector3D(xc',xc' + 0.5*element.n',[0 0 1]);
    view(-34,31);
end
function plotTriangle(i,p1,p2,p3)
line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'Color','black');
line([p2(1),p3(1)],[p2(2),p3(2)],[p2(3),p3(3)],'Color','black');
line([p3(1),p1(1)],[p3(2),p1(2)],[p3(3),p1(3)],'Color','black');
xc = (p1 + p2 + p3)/3;
%text(xc(1),xc(2),num2str(i));
end
function ht = vector3D(xi,xf,c)
x = xf - xi;
X = [xi(1,1), xf(1,1)];
Y = [xi(2,1), xf(2,1)];
Z = [xi(3,1), xf(3,1)];
hl = line( X, Y,Z, 'Color', c,'Linewidth',2);
hp = piramide(xi,xf,[0,0,1]);
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






