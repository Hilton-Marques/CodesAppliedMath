% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;

%Geometry Brannan
%Differential geometry William burke

A = [[5,4];[4,5]];
[V,D] = eig(A);
n = 100;
theta = linspace(0,2*pi,n);
c = [cos(theta);sin(theta)];
figure
hold on
axis equal
axis([-10,10,-10,10]);
e = A*c;
plot(e(1,:),e(2,:));
quiver(0,0,D(1,1)*V(1,1),D(1,1)*V(1,2));
quiver(0,0,D(2,2)*V(2,1),D(2,2)*V(2,2));
%p = [8;0];
p = e(:,1);
fac = diag(sqrt(2));
fac = 1;
e_trans = V*fac*D*V'*c + p;
plot(e_trans(1,:),e_trans(2,:));

[p1,p2] = circleIntersection(A\p);
pp = A*p1;
pp2 = A*p2;
plot(p(1),p(2),'o','markersize',8);
plot(pp(1),pp(2),'o')
plot(pp2(1),pp2(2),'o')
v = polar(p);
l = norm(v);
k = dot(p,v)/dot(v,v);
drawLine(v);
[perp1,perp2] = ortho(A,p);
f_ortho(p,perp1)

plot(perp1(1),perp1(2),'o');
plot(perp2(1),perp2(2),'o');

keyboard
function out = f(x)
A = [[5,4];[4,5]];
A_inv = inv(A);
out = dot(x,(A_inv'*A_inv)*x) - 1;
end
function out = f_ortho(u,v)
A = [[5,4];[4,5]];
A_inv = inv(A);
out = dot(u,(A_inv'*A_inv)*v);
end
function out = polar(x)
A = [[5,4];[4,5]];
A_inv = inv(A);
v = (A_inv'*A_inv)*x;
out = v/(dot(v,v));
end
function drawLine(v)
ortho = [-v(2);v(1)];
fac = 20;
pi = v + fac*ortho;
pj = v - fac*ortho;
line([pi(1),pj(1)],[pi(2),pj(2)]);
pi =  + fac*ortho;
pj =  - fac*ortho;
line([pi(1),pj(1)],[pi(2),pj(2)]);
end
function [p1,p2] = circleIntersection(x2)
x1 = [0;0];
r1 = 1;
r2 = sqrt(dot((x2 - x1),(x2 - x1)));
% figure
% hold on
%see wikipedia article on radical axis
dir = x2 - x1;
d = norm(dir);
n = dir/d;
a = (d^2 + r1^2 - r2^2) / (2*d);
b = sqrt(r1^2 - a^2);
n_ortho_1 = [-n(2); n(1)];
p1 = x1 + a*n + b*n_ortho_1;
p2 = x1 + a*n - b*n_ortho_1;
end
function [res1,res2] = ortho(A,p)
pi = A\p;
pj1 = [-pi(2);pi(1)];
pj2 = -[-pi(2);pi(1)];
res1 = A*pj1; multil
res2 = A*pj2; 
end