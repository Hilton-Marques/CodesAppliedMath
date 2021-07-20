clc;
clear;
close all;
%points = makeRand(1000);
points = xlsread('drill','Planilha1');

hold on
view(30,30);

[k1,av1] = convhull(points);

trisurf(k1,points(:,1),points(:,2),points(:,3),'FaceColor','cyan');
axis equal;
keyboard;
function points = makeRand(n)
points = zeros(n,3);
for i = 1:n
    points(i,:) = [rand,rand,rand];
end
end