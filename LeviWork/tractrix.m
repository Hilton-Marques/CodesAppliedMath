clc;
clear;
close all;

figure
hold on
dy = 0.01;
t = [-1;0];
R = norm(t);
ksi = (dy^2/R);
p0 = [1;0];
pts = [];
y = [0;1];
p_old = p0;
t_old = t;
for i = 0:100
    sin_teta = abs(dot(t,y));
    t_ortho = [t(2);-t(1)];
    t_ortho = t_ortho/norm(t_ortho);
    t_old = t;
    t = t + t_ortho*dy;
    p_old = p0;
    pps = [p_old,i*y*dy,(i+1)*y*dy];
    pps(3,:) = 0;
    trisurf([1,2,3],pps(1,:),pps(2,:),pps(3,:),'FaceColor','red');
    pps = [[1;0],[1;0]+t_old,[1;0]+t];
    pps(3,:) = 0;
    trisurf([1,2,3],pps(1,:),pps(2,:),pps(3,:),'FaceColor','blue');
    p0 = p0 + sin_teta*dy*t/(norm(t));    
    dh = p0 - p_old;
    dx = dh(1);
    -dh(1)/norm(dh);
    p_old(1)/R;
    %plot(p0(1),p0(2),'--','markersize',5);
    pts(:,end+1) = p0;
    pause(0.1);
end
%fimplicit(@(x,y) (1/sec(y)) - sqrt(1-y^2) - x)
plot(pts(1,:), pts(2,:))