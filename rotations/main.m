% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;


figure 
axis off
hold on
axis equal
axis([-0.2,3,-3,1.5])
set(gcf,'color','w');

xo = [0;1];
xf = [cos(0.4*pi);sin(0.4*pi)];

n = 100;
quiver(0,0,xo(1),xo(2),'color','blue','linewidth',1);
for i = 1:n
   h = i/n;
   theta = getAngle(xo) + h*(getAngle(xf) - getAngle(xo));
   tr = h*[2;-2];
   R = rot(theta);
   xi = R*[1;0];   
   h = quiver(tr(1),tr(2),xi(1),xi(2),'color','red','linewidth',1);   
   pause(0.01);
   delete(h);
end
h = quiver(tr(1),tr(2),xi(1),xi(2),'color','red','linewidth',1);


function T = trans(t)
T = eye(3);
T(1,3) = t(1);
T(2,3) = t(2);
end
function R = rot(theta)
R = eye(2);
R(1,1) = cos(theta);
R(2,1) = sin(theta);
R(1,2) = -sin(theta);
R(2,2) = cos(theta);
end
function theta = getAngle(u)
theta = mod(atan2(u(2),u(1)),2*pi);
end