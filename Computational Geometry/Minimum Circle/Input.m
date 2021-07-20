%% Clear memory
clc;
clear;
close all;
%points = xlsread('Input','Planilha1');
%points = points(1:end,1:2);
points = makeRand(200);
tic
circle = struct('radius',0,'center',[]);
%% Find all permutations
n = 1:1:length(points);
comb3 = nchoosek(n,3);
comb2 = nchoosek(n,2);
%% Start
for i = 1:length(comb2)
    p1 = points(comb3(i,1),:);
    p2 = points(comb3(i,2),:);
    circleTemp = makeCircle(p1,p2);
    if isAValidCircle(points,circleTemp)
        if (circleTemp.radius <= circle.radius || isempty(circle.center))
            circle = circleTemp;
        end
    end
end
for i = 1:length(comb3)
    p1 = points(comb3(i,1),:);
    p2 = points(comb3(i,2),:);
    p3 = points(comb3(i,3),:);
    circleTemp = makeCircle(p1,p2,p3);
    %drawCircle(circleTemp)
    if isAValidCircle(points,circleTemp)
        if (circleTemp.radius <= circle.radius || isempty(circle.center))
            circle = circleTemp;
        end
    end
end
toc
%% Plot
% hold on
% plot(points(:,1),points(:,2),'o');
% draw(points,circle);
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
function out = isAValidCircle(points,circle)
for i = 1:length(points)
    r = points(i,:) - circle.center;
    if (vecnorm(r) > circle.radius)
        out = false;
        return;
    end
end
out = true;
end
function draw(points,circle)
plot(points(:,1),points(:,2),'o');
drawCircle(circle);
end
function drawCircle(circle)
t = linspace(0,2*pi,100);
points = circle.radius*[cos(t);sin(t)] + circle.center';
plot(points(1,:),points(2,:));
end
function points = makeRand(n)
points = zeros(n,2);
for i = 1:n
    points(i,:) = [rand,rand]; 
end
end