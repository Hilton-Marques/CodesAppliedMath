clc
clear all;
close all;


% draw fast multipole
%triangle field
p1 = [-3, -1.5, 0.0 ];
p2 = [-2.25, -0.75 , 0.0];
p3 = [-2.25, -1.125, 1.0];
tri = Triangle(p1,p2,p3);
yc =  [-2, -0.75, 2/3];
xs =  [ 2.25, -0.75, 0.0 ];
xf = p1;
u = xf - yc;
v = xs - yc;
centroide = tri.getCentroide();
%plot
figure;
hold on
view(30,30);
tri.show();
plot3(yc(1),yc(2),yc(3),'o');
plot3(xs(1),xs(2),xs(3),'o');
quiver3( yc(1), yc(2), yc(3), u(1), u(2), u(3) ,1.0, 'r');
quiver3( yc(1), yc(2), yc(3), v(1), v(2), v(3) ,1.0,'r');
text(p1(1),p1(2),p1(3), 'P1');
text(p2(1),p2(2),p2(3), 'P2');
text(p3(1),p3(2),p3(3), 'P3');
text(yc(1),yc(2),yc(3), 'Xc', 'HorizontalAlignment', 'right',...
    'VerticalAlignment', 'top');
text(xs(1),xs(2),xs(3), 'Xs', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');
text(centroide(1),centroide(2),centroide(3), 'T');
set(gca,'visible','off')
