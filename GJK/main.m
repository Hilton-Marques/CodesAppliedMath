% Clear workspace
clear
close(findall(0,'Type','figure'));
clc

v1 = [-0.5,-0.5];
v2 = [0.5,-0.5];
v3 = [0.5,0.5];
v4 = [-0.5,0.5];
square1 = square([v1;v2;v3;v4]);
square2 = square([v1;v2;v3;v4]);
square1.transform(1*pi/2.5,1.0*[0.5,-1.0]);
%square2.transform(0*pi/2.5,1.0*[1.0,-1.0]);
circle1 = circle(1);
circle1.transform(0,-0.4*[1.5,-0.5]);
p_m = MikSub(square1,circle1);
k = convhull(p_m(:,1),p_m(:,2));

figure
hold on
axis equal;
set(gca,'visible','off');
set(gcf,'color','white');

circle1.plot('red');
square1.plot('blue');
plot(p_m(k,1),p_m(k,2),'color','black')
plot(0,0,'o','color','yellow','MarkerFaceColor','yellow');
%exportgraphics(gca,'foto4.jpeg','Resolution',1000);

[bool,point] = GJK(square1,circle1);
plot(point(1),point(2),'*');
%exportgraphics(gca,'foto2.jpeg','Resolution',1000);
%plot(p_m(k,1),p_m(k,2))
%plot(p_m(k,1),p_m(k,2),'*');

%plotBdMik(square1,square2);
%[bool,point] = GJK(square1,square2);



%plot(p_m(:,1),p_m(:,2),'*')
%plotBdMik(square1,circle1);

function p_m = MikSub(convex1,convex2)
n1 = size(convex1.pts,1);
n2 = size(convex2.pts,1);
p_m = zeros(n1*n2,2);
count = 1;
for i = 1:n1
    pi = convex1.pts(i,:);
    for j = 1:n2
        pj = convex2.pts(j,:);
        p_m(count,:) = pi - pj;
        count = count + 1;
    end
end
end
function plotBdMik(convex1,convex2)
n = convex1.n + convex2.n + 1;
teta = linspace(0,2*pi,n);
d = [cos(teta)',sin(teta)'];
bd = zeros(n,2);
for i = 1:n
    di = d(i,:);
    vA = convex1.suportfunction(di');
    vB = convex2.suportfunction(-di');
    bd(i,:) = vA - vB;
end
bd = unique(bd,'rows');
plot(bd(:,1),bd(:,2),'*');
%plot(bd(:,1),bd(:,2));

end
function [bool,p] = GJK(convex1,convex2)
%plot(0,0,'o','color','yellow');
triA = triangle();
triB = triangle();
tri = triangle();
d0 = (convex2.centroide - convex1.centroide);
d0 = d0/norm(d0);
[vA,vB] = support(convex1,convex2,d0);
A = vA - vB;
%plot(vA(1),vA(2),'*');
%plot(vB(1),vB(2),'*');
plot(A(1),A(2),'*','color','black');
exportgraphics(gca,'foto6.jpeg','Resolution',1000);
tri.append(A);
triA.append(vA);
triB.append(vB);
d = -A;
d = d/norm(d);
while true
    [vA,vB] = support(convex1,convex2,d);
    P = vA - vB;
    plot(P(1),P(2),'*','color','black');
    exportgraphics(gca,'foto8.jpeg','Resolution',1000);
    if dot(P,d) < 0
        bool = false;
        break
    end
    tri.append(P);
    triA.append(vA);
    triB.append(vB);
    [flag,d] = handleSimplex(tri,triA,triB);
    exportgraphics(gca,'foto9.jpeg','Resolution',1000);
    if flag
        bool = true;
        break;
    end
    
end
lam = tri.getOriginBary();
tri.plot('green');
exportgraphics(gca,'foto10.jpeg','Resolution',1000);
triA.plot('red');
triB.plot('blue');
%exportgraphics(gca,'foto11.jpeg','Resolution',1000);

p = triA.plotBaryPoint(lam);
exportgraphics(gca,'foto7.jpeg','Resolution',1000);

p = triB.plotBaryPoint(lam);
end
function [vA,vB] = support(convex1,convex2,di)
    vA = convex1.suportfunction(di');
    vB = convex2.suportfunction(-di');
    v =  vA - vB;
    d_unit = 0.4*di/norm(di);
    quiver(convex1.centroide(1),convex1.centroide(2),d_unit(1),d_unit(2),'color','magenta');
    quiver(convex2.centroide(1),convex2.centroide(2),-d_unit(1),-d_unit(2),'color','magenta');
    plot(vA(1),vA(2),'o');
    plot(vB(1),vB(2),'o');
end

function [bool, d] = handleSimplex(tri,triA,triB)
if tri.n == 2
    [bool,d] = tri.GetPerpDir2Line();
    return
end
[bool, d] = tri.ContainOrigin(triA,triB);
end
function MerryGoRound()
n = 5;
t = linspace(0,2*pi,n);
for i = 1:n
    square1.plot('blue');
    square2 = square([v1;v2;v3;v4]);
    square2.transform(t(i),1.0*[1.0,-1.0]);
    square2.plot('red');
    p_m = MikSub(square1,square2);
    k = convhull(p_m(:,1),p_m(:,2));
    plot(p_m(k,1),p_m(k,2));
    [bool,point] = GJK(square1,square2);
    plot(point(1),point(2),'*');
    pause(1);
    clf
    figure
    hold on
end
end

