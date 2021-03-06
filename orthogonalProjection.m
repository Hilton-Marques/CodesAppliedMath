clc;
clear;
close all;
[X,Y,Z] = sphere();
hold on
view(30,30)
axis equal;
set(gca,'visible','off');
set(gcf,'color','white');
surf(X,Y,Z,'FaceAlpha',0.5)

drawPlan()
c1 = drawCurve([0,0],[pi/2,0]);
c2 = drawCurve([pi/2,0],[pi/2,pi/2]);
pause(0.5);
projectCurve(c1);
projectCurve(c2);

axis equal

function drawPlan()
A = 2.5;
color = [1,0,0];
xo = [0,0,-2];
x = [1,0,0];
y = [0,1,0];
xp1 = xo + A*x;
yp1 = xo + A*y;
xp2 = xo - A*x;
yp2 = xo - A*y;
h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.5);
end
function c = drawCurve(A,B)
teta1 = A(1);
fi1 = A(2);
teta2 = B(1);
fi2 = B(2);
tetas = linspace(teta1,teta2,20);
fis = linspace(fi1,fi2,20);
x = sin(tetas) .* cos(fis);
y = sin(tetas) .* sin(fis);
z = cos(tetas);
c = [x',y',z'];
plot3(x,y,z,'color','blue','linewidth',5);
end
function projectCurve(c1)
x = c1(:,1);
y = c1(:,2);
z = c1(:,3);
pts_projected = [];
for i = 1:size(x,1)
    x_i = x(i,1);
    y_i = y(i,1);
    z_i = z(i,1);
    ray_begin = [x_i;y_i;2];
    ray_end = [x_i;y_i;-2];
    line([ray_begin(1), ray_end(1)], [ray_begin(2), ray_end(2)], [ray_begin(3), ray_end(3)],'color','yellow');
    pts_projected(end+1,1:3) = ray_end;
    pause(0.5)
    plot3(pts_projected(:,1),pts_projected(:,2),pts_projected(:,3),'color',[0.5,0.5,0.5],'linewidth',5);
end
end