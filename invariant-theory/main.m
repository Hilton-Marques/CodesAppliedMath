clc;
clear;
close all;
rng('default');
A = [[-2,1];[-3,2]];
A = [[0,-1];[1,0]];

figure
hold on

n = 100;
x = linspace(-1,1, n);
[X,Y] = meshgrid(x);
[U,V] = chi(X,Y);
quiver(X,Y,U,V);
Z = phi(X,Y);
contour(X,Y,Z,200);
x = X(:);
y = Y(:);

x = rand(10,1)*2 -1;
y = rand(10,1)*2 -1;
colors = rand(10,3);
scatter(x,y,30,colors,'filled');
for i = 1:4
    v = A^i * [x';y'];
    scatter(v(1,:),v(2,:),30,colors,'filled');
end
%surf(X,Y,Z)
exportgraphics(gca,'rotation-example.jpg','Resolution',300)
function [U,V] = chi(X,Y)
% U = -2*X + Y;
% V = -3*X + 2*Y;
U = -Y;
V = X;
%U = -2 .* Y .* X .^2;
%V = 2 .* X .* Y .^2;
end
function Z = phi(X,Y)
%Z = X .* (Y - 2*X);
Z = X.^2 .* Y.^2;
end