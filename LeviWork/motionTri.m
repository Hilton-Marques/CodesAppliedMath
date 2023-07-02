clc;
clear;
close all;

figure
hold on
axis equal
view(0,90);


p0 = [0;0;0];
p1 = [1;0;0];
p2 = [1;1;0];
pts_2 = [p0,p1,p2];


p0 = [0;0;0];
p1 = [1;0;0];
p2 = [0;1;0];
pts = [p0,p1];
setBB(pts');
quiver(0,0,p1(1),p1(2),0);
angle = 2*pi;
n = 10;
t = linspace(0,angle,n);
trans = [1;1];
p_old = p2;
p_lin = p1;
for i = 1:n
    teta = t(i);
    R = Rot(teta);
    %p_new = R*p1;
    ptsi = R*pts_2;

    ortho = [-p_old(2);p_old(1);0];
    ortho = ortho/norm(ortho);
    p_new = p_old + angle/(n)*ortho;
    p_dir = (p_old/norm(p_old)) * angle/(n);
    p_lin = p_lin + p_dir;
    DrawTri([p_lin,p_lin + p_old, p_lin + p_new]');
    p_old = p_new;
    
    p_i_n = ptsi(:,2);
    p_j_n = ptsi(:,3);
    f = DrawTri([p0,p_i_n,p_j_n]');
    pause(0.01);
    delete(f);
end
DrawTri([p0,p_i_n,p_j_n]');
%moveTri(angle);
% p_f = Rot(angle) * p1;
% quiver(0,0,p_f(1),p_f(2),0);



function moveTri(angle)
p0 = [0;0;0];
p1 = [1;0;0];
p2 = [1;1;0];
pts = [p0,p1,p2];
% figure
% hold on
% axis equal
% view(0,90);
% setBB(pts);

n = 120;
t = linspace(0,angle,n);
for i = 1:n
    teta = t(i);
    R = Rot(teta);
    ptsi = R*pts;
    p_i_n = ptsi(:,2);
    p_j_n = ptsi(:,3);
    f = DrawTri([p0,p_i_n,p_j_n]');
    pause(0.01);
    delete(f);
end
end

function f = DrawTri(pts)
f = trisurf([1,2,3],pts(:,1),pts(:,2),pts(:,3),'FaceAlpha',1.0,'EdgeAlpha',1.0);
end
function setBB(pts)
margin = 1.5;
x_min = min(pts(:,1)) - margin;
x_max = max(pts(:,1)) + margin;
y_min = min(pts(:,2)) - margin;
y_max = max(pts(:,2)) + margin;
z_min = min(pts(:,3)) - margin;
z_max = min(pts(:,3)) + margin;
bb = [x_min,x_max,y_min,y_max,z_min, z_max];
axis(bb);
end
function R = Rot(teta)
R = [[cos(teta),-sin(teta),0];[sin(teta),cos(teta),0];[0,0,1]];
end
function T = rigidMotion(teta,t)
R = Rot(teta);
T = Trans(t);
T = T * R;
end
function T = Trans(t)
T = eye(3);
T(1,3) = t(1);
T(2,3) = t(2);
end