clc;
clear;
close all;

hold on
view(30,30)
axis([0,10,0,10,0,5]);
%set(gca,'visible','off');
set(gcf,'color','white');
drawPlan()
x = linspace(0,2.5,10);
y = x.^2 ;
pts = [x;y];
%plot3(x,y,ones(10,1));
%plotRay(pts);
plotMarch();
function drawPlan()
A = 25;
color = [1,0,0];
xo = [0,0,1];
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
function plotRay(pts)
n = size(pts,2);
pts_z = [pts;ones(1,n)];
alpha = linspace(0,10,10);
for i = 1:10
    p0 = [0;0;0];
    projected_pts = alpha(i)*pts_z;
    for j = 1:n
      pj = projected_pts(:,j);
      line([p0(1),pj(1)],[p0(2),pj(2)],[p0(3),pj(3)],'color','yellow');
      plot3(projected_pts(1,1:j),projected_pts(2,1:j),projected_pts(3,1:j));
    end
    pause(1);
end

end
function plotMarch()
    n = 100;
    pi = [1;1];
    alpha = linspace(0.5,n,100);
    r = [0;1];
    pj = pi + alpha.*r;
    pj_z = [pj;ones(1,n)];
    for i = 1:100
        alpha_i = alpha(i);        
        pi_0 = pj_z(:,1);
        pi = pj_z(:,i);
        p0 = [0;0;0];
        projected_p0 = (1/alpha_i)*pj_z(:,1);
        projected_pi = (1/alpha_i)*pj_z(:,i);
        fac = 5;
        line([p0(1),fac*projected_p0(1)],[p0(2),fac*projected_p0(2)],[p0(3),fac*projected_p0(3)],'color','yellow');
        line([p0(1),fac*projected_pi(1)],[p0(2),fac*projected_pi(2)],[p0(3),fac*projected_pi(3)],'color','yellow');
        
        line([projected_p0(1),projected_pi(1)],[projected_p0(2),projected_pi(2)],[projected_p0(3),projected_pi(3)],'color','red');
        
        line([pi_0(1),pi(1)],[pi_0(2),pi(2)],[pi_0(3),pi(3)],'color','black');
        
        pause(1);
    end
end