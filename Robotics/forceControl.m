clc;
clear all;
close all;

hold on
set(gca,'visible','off');
view(43,24);
axis equal
xp1 = [0,0,0];
yp1 = [3,0,0];
xp2 = [3,1,0];
yp2 = [0,1,0];
color = [0.8,0.8,0.8];

h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.5);

xp1 = [0,0.5,0];
yp1 = [3,0.5,0];
xp2 = [3,0.8,0];
yp2 = [0,0.8,0];
color = [1,0,0];

h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.5);
p = [0,0.5,0.1];
drawXY(p,[1,0,0],[0,0,1]);
drawZ(p,[0,0,1]);
x1 = [0,0.5,0];
x2 = [0.4,0.5,1];
h1 = line([x1(1),x2(1)],[x1(2),x2(2)],[x1(3),x2(3)],'Color',[0,0.5,0.3],...
    'Linewidth',10);
function drawPlan(n,h,w,color,xo)
if nargin == 1
    h = 1;
    w = 1;
    xo = [0;0;0];
    color = [1,0,0];
elseif nargin == 3
    xo = [0;0;0];
    color = [1,0,0];
elseif nargin == 4
    xo = [0;0;0];
end
[x,y] = findTriedro(n);
xp1 = xo + h*x;
yp1 = xo + w*y;
xp2 = xo - h*x;
yp2 = xo - w*y;
h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.5);
end
function [x,y] = findTriedro(z)
z = z/norm(z);
uTemp = [1;1;1];
uTemp = uTemp/norm(uTemp);
if (dot(z,uTemp) == 1)
    uTemp = uTemp([2,1,3]);
end
x = cross(uTemp,z)/norm(cross(uTemp,z));
y = cross(z,x)/(norm(cross(z,x)));
end
function drawXY(p,x,z)
y = cross(z,x);
posTextX = p + x;
posTextY = p + y;
text(posTextX(1),posTextX(2),posTextX(3),strcat('x'));
text(posTextY(1),posTextY(2),posTextY(3),strcat('y'));
quiver3(p(1),p(2),p(3),...
    x(1),x(2),x(3),'Color','red');
quiver3(p(1),p(2),p(3),...
    y(1),y(2),y(3),'Color','green');

end
function drawZ(p,z)
quiver3(p(1),p(2),p(3),...
    z(1),z(2),z(3),'Color','blue');
posTextZ = p + z;
text(posTextZ(1),posTextZ(2),posTextZ(3),strcat('z'),...
    'HorizontalAlignment','right');
end