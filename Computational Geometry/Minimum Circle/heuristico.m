clc;
clear;
close all;
%points = xlsread('Input','Planilha1');
points = makeRand(200);
%% Find further points
tic
[~,index1] = min(points(:,1));
[~,index2] = min(points(:,2));
[~,index3] = max(points(:,1));
[~,index4] = max(points(:,2));
p1 = points(index1,:);
p2 = points(index2,:);
p3 = points(index3,:);
p4 = points(index4,:);
p = [p1;p2;p3;p4];
d = 0;
for i = 1:4
    pi = p(i,:);
    for j = i+1:4
        pj = p(j,:);
        dTemp = vecnorm(pi - pj);
        if dTemp > d
            d = dTemp;
            fp = [pi;pj];
        end
    end
end
circle = makeCircle(fp(1,:),fp(2,:));
%% Update circle
for k = 1:length(points)
    pk = points(k,:);
    r = (pk - circle.center);
    d = vecnorm(r);
    if ( d > circle.radius)
        circle.center = circle.center + (0.5/d)*(d - circle.radius)*r;
        circle.radius = (d + circle.radius)/2;
    end
end
toc
%% Plot
% hold on
% plot(points(:,1),points(:,2),'o');
% drawCircle(circle);

function circle = makeCircle(p1,p2)
circle.center = (p1 + p2)/2;
circle.radius = vecnorm(p1 - p2)/2;
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