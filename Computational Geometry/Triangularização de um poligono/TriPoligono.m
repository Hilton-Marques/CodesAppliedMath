clc;
clear;
close all;
%% Solução de Skiena & Revilla para tecer um polígono
points = xlsread('Input','Planilha1');
N = size(points,1);
%% StartTriangularization
obj = TriPolygon(points,N);
obj.startTri();
%% Export and Plot
fileID = fopen('indexTri.txt','w');
hold on
for i = 1:N
p = [points(i,1),points(i,2)];
plot(p(1),p(2),'o','Color','Red');
text(p(1),p(2),num2str(i))
end
for i = 1:length(obj.triIndex)
    j = obj.triIndex(i,1);
    k = obj.triIndex(i,2);
    l = obj.triIndex(i,3);
    plotTriangle(i,points(j,:),points(k,:),points(l,:));
    fprintf(fileID,'%5d %5d %5d\n',[j,k,l]);
end
function plotTriangle(i,p1,p2,p3)
    line([p1(1),p2(1)],[p1(2),p2(2)],'Color','blue');
    line([p2(1),p3(1)],[p2(2),p3(2)],'Color','blue');
    line([p3(1),p1(1)],[p3(2),p1(2)],'Color','blue');
    xc = (p1 + p2 + p3)/3;
    text(xc(1),xc(2),num2str(i));
end
