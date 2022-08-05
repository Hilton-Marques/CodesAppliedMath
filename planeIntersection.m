% Clear workspace
clear
close(findall(0,'Type','figure'));
clc


intersect3D([1,1,6],[3,2,5])
keyboard
l1 = 1*[1,1];
l2 = [3,2];
figure
hold on
plot(0,0,'o','markersize',10)
drawLine(l1);
drawLine(l2)
%intersect(l1,l2)

function drawLine(l)
lam = 22;
a = l;
ortho1 = [-l(2),l(1)];
p1 = a + lam*ortho1/norm(ortho1);
p2 = a - lam*ortho1/norm(ortho1);
quiver(0,0,l(1),l(2));
line([p1(1),p2(1)],[p1(2),p2(2)]);
end
function p = intersect(n,m)
c = dot(n,n);
d = dot(m,n);
m_perp = [-m(2),m(1)];
e = dot(m_perp,n);
lam2 = (c-d)/e;
p = m + lam2*m_perp;
plot(p(1),p(2),'o','markersize',10)
end
function p = intersect3D(n,m)
figure 
hold on
view(30,30);
quiver3(0,0,0,n(1),n(2),n(3));
quiver3(0,0,0,m(1),m(2),m(3));
z = cross(m,n);
m_perp = cross(z,m);
c = dot(n,n);
d = dot(m,n);
e = dot(m_perp,n);
lam2 = (c-d)/e;
p = m + lam2*m_perp;
quiver3(n(1),n(2),n(3),p(1)-n(1),p(2)-n(2),p(3)-n(3));
quiver3(m(1),m(2),m(3),p(1)-m(1),p(2)-m(2),p(3)-m(3));
lam = 5;
p1 = p + lam*z;
p2 = p - lam*z;
line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)]);
check1 = dot(p-n,n);
check2 = dot(p-m,m);
end