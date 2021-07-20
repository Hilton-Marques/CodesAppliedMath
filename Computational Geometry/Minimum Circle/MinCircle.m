clc;
clear;
close all;
%points = xlsread('Input','Planilha1');
points = makeRand(200);
tic
%% Random permutation
n = length(points);
index = randperm(n);
points = points(index,:);
%% Initialize first circle
circle = makeCircle(points(1,:),points(2,:));
%% Start
for i = 3:n
    pi = points(i,:);
    if (~PointBelongToCircle(circle,pi))
        circle = MinCircleWithPoint(points(1:i-1,:),pi);
    end
end
toc
%% Plot 
% hold on
% plot(points(:,1),points(:,2),'o');
% drawCircle(circle)
% % circle2.radius = 200.56;
% % circle2.center = [409.95,346.51];
% % drawCircle(circle2);
% legend({'Pontos','Círculo Mínimo','Círculo Heurístico'})
function circle = makeCircle(p1,p2,p3)
if nargin == 2
    circle.center = (p1 + p2)/2;
    circle.radius = vecnorm(p1 - p2)/2;
else
    l1 = p2 - p1;
    l2 = p3 - p1;
    a = l1;
    b = l2;
    area = cross([a,0],[b,0]);
    if (area == 0)
        a(1) = min([p1(1,1),p2(1,1),p3(1,1)]);
        a(2) = min([p1(1,2),p2(1,2),p3(1,2)]);
        b(1) = max([p1(1,1),p2(1,1),p3(1,1)]);
        b(2) = max([p1(1,2),p2(1,2),p3(1,2)]);
        circle.center = (b + a)/2;
        circle.radius = vecnorm(b - a)/2;
    else
        if (area < 0)
            a = l2;
            b = l1;
        end
        ac = [-a(2),a(1)];
        bc = [-b(2),b(1)];
        rc = (dot(b,b)*ac - dot(a,a)*bc)/(2*dot(ac,b));
        circle.center = p1 + rc;
        circle.radius = vecnorm(rc);
    end
end
end
function bool = PointBelongToCircle(circle,pi)
r = (pi - circle.center);
bool = true;
if (vecnorm(r) > circle.radius)
    bool = false;
    return
end
end
function circle = MinCircleWithPoint(points,q)
circle = makeCircle(points(1,:),q);
for j = 2:size(points,1)
    pj = points(j,:);
    if (~PointBelongToCircle(circle,pj))
        circle = MinCircleWith2Points(points(1:j-1,:),pj,q);
    end
end
end
function circle = MinCircleWith2Points(points,q1,q2)
circle = makeCircle(q1,q2);
for k = 1:size(points,1)
    pk = points(k,:);
    if (~PointBelongToCircle(circle,pk))
        circle = makeCircle(pk,q1,q2);
    end
end
end
function drawCircle(circle,p)
t = linspace(0,2*pi,100);
points = circle.radius*[cos(t);sin(t)] + circle.center';
plot(points(1,:),points(2,:));
if nargin > 1
    plot(p(:,1),p(:,2),'o','Color','red');
end
end
function points = makeRand(n)
points = zeros(n,2);
for i = 1:n
    points(i,:) = [rand,rand]; 
end
end