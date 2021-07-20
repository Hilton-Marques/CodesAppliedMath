clc;
clear;
close all;
%% Solução Bowyer – Watson 
points = xlsread('Input','Planilha4');
%points = [rand([5,1]),rand([5,1])];
points = [0.9734, 0.9224;
0.6441,0.1656;
0.9813,0.6273;
0.1804,0.6160;
0.6765,0.6146;
0.3936,0.8655;
0.5042,	0.0374;
0.1221,	0.0343;
0.1528,	0.6141;
0.1175,	0.7539];
hold on
for i = 1:length(points)
    p = points(i,:);
    plot(p(1),p(2),'o','Color','blue');
end
hold off
N = size(points,1);
%[~,I] = sort(points(:,1));
%points = points(I,:);
DT = delaunay(points(:,1),points(:,2));
figure
hold on
triplot(DT,points(:,1),points(:,2));
hold off
%% Find super Triangle
xmax = [max(points(:,1)),max(points(:,2))];
xmin = [min(points(:,1)),min(points(:,2))];
xMed = (xmax + xmin)/2;
M1 = (xmax(1,1) - xmin(1,1))/2;
M2 = (xmax(1,2) - xmin(1,2))/2;
M = max(M1,M2);
p1 = xMed - 3*[M,M];
p2 = xMed + [0,3*M];
p3 = xMed + [3*M,0];
neigh = [Triangle,Triangle,Triangle];
triangleRoot = Triangle(p1,p2,p3,neigh); 
for i = 1:3
    triangleRoot.edges(i).T = triangleRoot;
end
DT = DelaunayTriangulation(triangleRoot);
%% StartTriangularization
for i = 1:N
    pi = points(i,:);
    DT.addPoint(pi);
end
%% Delete edges from the superTriangle
figure
hold on
triangles = Triangle.empty;
pointsSuper = [p1;p2;p3];
for i = 1:length(DT.triangles)
    Ti = DT.triangles(i);
    Ti.plot([0,0,1]);
    pMed = (Ti.points(1,:)+Ti.points(2,:)+Ti.points(3,:))/3;
    text(pMed(1),pMed(2),num2str(i));
    index = ~ismember(pointsSuper,Ti.points,'rows'); 
    if (sum(index) == 3)
        triangles(end + 1) = Ti;
    end
end
hold off
T = length(triangles);
figure
hold on
for i = 1:T
    Ti = triangles(i);
    Ti.plot([0,0,1]);
end
