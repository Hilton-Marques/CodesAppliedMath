clc;
clear;
close all;

addpath("../../my_libs/cgeom/");

p0 = [0,0,0];
p1 = [1,0,0];
p2 = [0,1,0];
p3 = [0,0,1];
t1 = [p0;p1;p2];
t2 = [p1;p0;p3];
% t2 = t2 + n2*3 + n1*1;
% t1 = t1 + n1*1 + n2*3;

t1 = [0.5371    8.5351   -0.0073;...
    0.5356    8.5423   -0.0016;...
    0.5403    8.5418   -0.0069];

t2 =[0.5371    8.5351   -0.0073;...
    0.5307    8.5348   -0.0015;...
    0.5356    8.5423   -0.0016];

n1 = cgeom.normal(t1);
n2 = cgeom.normal(t2);
n1 = n1/norm(n1);
n2 = n2/norm(n2);


c = dot(n1,t1(1,:));
d = dot(n2,t2(1,:));
[o,d] = PlanePlaneIntersection(c*n1,n2*d);

sin_theta = norm(cross(n2,n1));
cos_theta = dot(n1,n2);

theta = atan2(sin_theta, cos_theta);

R = expMap(d*theta);

p_teste = (t2(2,:) - o)*R' + o;
%plots
hold on
set(gcf,'color','white');
            lighting gouraud;
            camlight('headlight');
            view(-9,48);
trisurf([1,2,3],t1(:,1),t1(:,2),t1(:,3),'FaceColor','red');
trisurf([1,2,3],t2(:,1),t2(:,2),t2(:,3),'FaceColor','blue');
view(30,30);

%plot3(p_teste(1),p_teste(2),p_teste(3),'o','markersize',5,'markerfacecolor','red');
close_pt = cgeom.closestPointToTriangle(t1, p_teste);
close_pt_2 = cgeom.closestPointToTriangle(t1, t2(2,:));

%plot3(close_pt(1),close_pt(2),close_pt(3),'o','markersize',5,'markerfacecolor','green');
%plot3(p_teste(1),p_teste(2),p_teste(3),'o','markersize',5,'markerfacecolor','red');

%plot3(t2(3,1),t2(3,2),t2(3,3),'o','markersize',5,'markerfacecolor','red');
check = dot(p_teste - o,n1);

%% check edges
s = [t1(1,:);t1(2,:)];
plot3(s(1,1),s(1,2),s(1,3),'o','markersize',5,'markerfacecolor','cyan');
plot3(s(2,1),s(2,2),s(2,3),'o','markersize',5,'markerfacecolor','magenta');
projected_pt = cgeom.project(t1, t2(2,:));
r = [projected_pt; t1(3,:)];
p = cgeom.segmentSegmentIntersection(r,s);
plot3(p(1),p(2),p(3),'o','markersize',5,'markerfacecolor','black');

keyboard;

function R = expMap(v)
theta = norm(v);
if (theta ~= 0)
    v = v/theta;
end
U = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];
R = eye(3) + sin(theta)*U + (1 - cos(theta))*U*U;
end

function [o,d] = PlanePlaneIntersection(n,m)
z = cross(m,n);
m_perp = cross(z,m);
c = dot(n,n);
d = dot(m,n);
e = dot(m_perp,n);
lam2 = (c-d)/e;
o = m + lam2*m_perp;
lam = 10000;
p1 = o + lam*z/norm(z);
p2 = o - lam*z/norm(z);
d = z/norm(z);
%line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'linewidth',3,'color','black');
check1 = dot(o-n,n);
check2 = dot(o-m,m);
end
