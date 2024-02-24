clear;
clc;
close all;

%Q = eye(2);
Q = [[5  ,4];[4, 5]];
det(Q)
% d = @(x,y) dot(x, Q*y);
%d = @(x,y) diag((y - x)' * Q * (y - x));
d = @(x,y) (y - x)' * Q * (y - x);

x = [0;0];


figure
hold on
% n = 100;
% xi = linspace(x(1),y(1),n);
% yi = linspace(x(2),y(2),n);
% r = 0.5;
% [X,Y] = meshgrid(xi,yi);

% pts = [X(:), Y(:)]';
% Z = cost(r, pts);
% Z = reshape(Z,[n,n]);
% Z = Z / max(Z(:));
% [m n] = min(Z(:));
% [x y] = ind2sub(size(Z),n);
% X(x,y)
% Y(x,y)
% 
% r = 18; % #levellines
% clf; hold on;
% imagesc(xi, yi, Z');
% contour(xi,yi,Z',linspace(0,1,r), 'k');
% % contour(m,s,f(m1,s1,M,S)',f(m1,s1,m1*ones(1,r),linspace(1e-2,s1,r)), 'k');
% plot(x(1),x(2), 'k.', 'MarkerSize', 20);
% colormap(parula(r-1));
% % colormap(parula(256));
% caxis([0 1]);
% axis equal; axis off;
% drawnow;
% 
% keyboard
options = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','quasi-newton','MaxIterations',1);
theta = linspace(0,2*pi,10);
for i = 1:10
   %y = [2;-2];
   y = [cos(theta(i)); sin(theta(i))];

cost = @(r, m) (1 - r) * d(x , m) + r * d(y , m);
for r = 0:0.1:1
    h = @(m) cost(r, m);
    [res, val] = fminunc(h,[1;1]);
    plot(res(1),res(2),'o');
    %delete(graph);
end
end
figure 
hold on
Q = eye(2);
d = @(x,y) (y - x)' * Q * (y - x);
cost = @(r, m) (1 - r) * d(x , m) + r * d(y , m);

for r = 0:0.01:1
    h = @(x) cost(r, x);
    [res, val] = fminunc(h,[1;1]);
    plot(res(1),res(2),'o');
    %delete(graph);
end