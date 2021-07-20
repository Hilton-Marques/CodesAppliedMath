clc;
clear all;
close all;
%% Input data
% %% MIT DirectDrive Arm
% nBarras = 5;
% nJuntas = 3;
% l1 = 2;
% l2 = 1;
% l3 = 3;
% l4 = 0.8;
% l5 = 3;
% link1 = [0;0;1];
% link2 = [1;0;0];
% link3 = [0;1;0];
% link4 = [0;0;1];
% link5 = [1;0;0];
% z0 = [0;0;1];
% z1 = [0;1;0];
% z2 = [0;0;1];
% L = [l1,l2,l3,l4,l5];
% members = [link1,link2,link3,link4,link5];
% joints = [z0,z1,z2];
% O = [0;0;0];
% positions = zeros(3,nBarras+1);
% indexJuntas = [1,3,4];
% flag_JuntasOff = true;
%%  Schilling Titan II
% nBarras = 6;
% nJuntas = 6;
% l1 = 1;
% l2 = 1;
% l3 = 3;
% l4 = 3;
% l5 = 1;
% l6 = 1.5;
% O = [0;0;0];
% link1 = [0;0;1];
% link2 = [0;1;0];
% link3 = [0;0;1];
% link4 = [0;1;0];
% link5 = [0;0;-1];
% link6 = [1;0;0];
% z0 = [0;0;1];
% z1 = [1;0;0];
% z2 = [1;0;0];
% z3 = [1;0;0];
% z4 = [0;1;0];
% z5 = [1;0;0];
% L = [l1,l2,l3,l4,l5,l6];
% members = [link1,link2,link3,link4,link5,link6];
% joints = [z0,z1,z2,z3,z4,z5];
% positions = zeros(3,nBarras+1);
% indexJuntas = [1,2,3,4,5,6];
% flag_JuntasOff = false;
%% Prova
nBarras = 6;
nJuntas = 6;
l1 = 1.1;
l2 = 0.9;
l3 = 1.3;
l4 = 0;
l5 = 0;
l6 = 1.5;
O = [0;0;0];
link1 = [0;0;1];
link2 = [0;0;1];
link3 = [0;1;0];
link4 = [0;1;0];
link5 = [1;0;0];
link6 = [1;0;0];
z0 = [0;0;1];
z1 = [0;0;1];
z2 = [0;1;0];
z3 = [0;1;0];
z4 = [0;0;1];
z5 = [1;0;0];
L = [l1,l2,l3,l4,l5,l6];
members = [link1,link2,link3,link4,link5,link6];
joints = [z0,z1,z2,z3,z4,z5];
positions = zeros(3,nBarras+1);
indexJuntas = [1,2,3,4,5,6];
flag_JuntasOff = [true,false,false,false,false,false];
juntasType = [1,0,0,1,1,1]; % 1 for rotational 0 for prismatic
%% Plot
hold on
set(gca,'visible','off');
view(43,24);
axis equal
% drawBase
drawPlan(link1,1.5);
x = [1;0;0];
z = [0;0;1];
if (L(1) ~= 0)
    drawZ(O,z,0);
    drawXY(O,x,z,0);
end
%% Draw Members
p0 = O;
for i = 1:nBarras
    positions(:,i) = p0;
    p1 = p0 + L(i)*members(:,i);
    drawMember(p0,p1);
    p0 = p1;
end
positions(:,end) = p1;
%% Draw Joints
for i = 1:nJuntas
    index = indexJuntas(i);
    LJoint = L(index);
    if (LJoint == 0)
        continue
    end
    if (flag_JuntasOff(i) == false)
        LJoint = 0;
    end
    z = joints(:,i);
    pJoint = positions(:,index);
    if juntasType(i)
        drawJoint(pJoint,z,0.5*LJoint);
    else
        drawPrism(pJoint,z,0.5*LJoint);
    end
end
%% Draw Coordinate System
for i = 2:nJuntas
    index = indexJuntas(i);
    pOrigin = positions(:,index)
    z = joints(:,i)
    zBefore = joints(:,i-1);
    drawZ(pOrigin,z,i-1);
    x = -cross(z,zBefore);
    if x == 0
    x = [1;0;0];
    end
    x
    drawXY(pOrigin,x,z,i-1);
end
%% EndEffect
zf = joints(:,nJuntas);
pGarra = positions(:,end)
drawPlan(members(:,end),1,pGarra,[.8,.8,.8]);
drawZ(pGarra,zf,nJuntas);
x = x;
drawXY(pGarra,x,zf,nJuntas);



%% Functions
function drawPlan(n,A,xo,color)
if nargin == 1
    A = 1;
    xo = [0;0;0];
    color = [1,0,0];
elseif nargin == 2
    xo = [0;0;0];
    color = [1,0,0];
end
[x,y] = findTriedro(n);
xp1 = xo + A*x;
yp1 = xo + A*y;
xp2 = xo - A*x;
yp2 = xo - A*y;
h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.5);
end
function drawMember(x0,x1)
line([x0(1),x1(1)],[x0(2),x1(2)],[x0(3),x1(3)],'Color',[0.5,0.5,0.5],...
    'Linewidth',20);
end
function drawJoint(x0,axis,offset)
if nargin == 2
    offset = 0;
end
x0 = x0 + offset*axis;
%axis = axis/norm(axis);
r = 0.2;
h = 0.3;
axisRot = cross([0,0,1],axis);
if (norm(axisRot) == 0)
    [X,Y,Z] = cylinder(r);
    Z = Z - 0.5;
    h2 = surf(X,Y,h*Z,'EdgeColor','none','facealpha',0.5);
    h2.XData = h2.XData + x0(1);
    h2.YData = h2.YData + x0(2);
    h2.ZData = h2.ZData + x0(3);
    return
end
[X,Y,Z] = cylinder(r);
Z = Z - 0.5;
h2 = surf(X,Y,h*Z,'EdgeColor','none','facealpha',0.5);
rotate(h2,axisRot,90,[0,0,0]);
h2.XData = h2.XData + x0(1);
h2.YData = h2.YData + x0(2);
h2.ZData = h2.ZData + x0(3);
end
function drawPrism(x0,axis,offset)
if nargin == 2
    offset = 0;
end
[x,y] = findTriedro(axis);
L = 0.6;
h = 0.01;
x1 = x0 + offset*axis - h*x;
x2 = x1 + L*axis;
x3 = x0 + offset*axis + h*x;
x4 = x3 + L*axis;
h1 = line([x1(1),x2(1)],[x1(2),x2(2)],[x1(3),x2(3)],'Color',[0,0.5,0.3],...
    'Linewidth',25);
h2 = line([x3(1),x4(1)],[x3(2),x4(2)],[x3(3),x4(3)],'Color',[0,0.5,0.3],...
    'Linewidth',25);


end

function [x,y] = findTriedro(z)
z = z/norm(z);
uTemp = [1;1;1];
uTemp = uTemp/norm(uTemp);
if (dot(z,uTemp) == 1)
    uTemp = uTemp([2,1,3]);
end
x = cross(uTemp,z);
y = cross(z,x);
end
function drawZ(p,z,index)
quiver3(p(1),p(2),p(3),...
    z(1),z(2),z(3),'Color','blue');
posTextZ = p + z;
text(posTextZ(1),posTextZ(2),posTextZ(3),strcat('z',num2str(index)),...
    'HorizontalAlignment','right');
end
function drawXY(p,x,z,index)
y = cross(z,x);
posTextX = p + x;
posTextY = p + y;
text(posTextX(1),posTextX(2),posTextX(3),strcat('x',num2str(index)));
text(posTextY(1),posTextY(2),posTextY(3),strcat('y',num2str(index)));
quiver3(p(1),p(2),p(3),...
    x(1),x(2),x(3),'Color','red');
quiver3(p(1),p(2),p(3),...
    y(1),y(2),y(3),'Color','green');
text(p(1)+0.2*x(1),p(2),p(3),strcat('O',num2str(index)),'HorizontalAlignment','right');

end