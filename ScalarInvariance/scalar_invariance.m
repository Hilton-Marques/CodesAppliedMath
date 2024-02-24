clc;
clear;
close all;



x = linspace(-10,10,100);
[X,Y] = meshgrid(x);
Z = f(X,Y);
[C,h] = contour(X,Y,Z,10);
plotl([1;0]);
plotl([-2;1]);
axis([-10,10,-10,10])
%clabel(C,h)
function Z = f(X,Y)
%Z = (X .* Y) ./ (X.^2 + Y.^2);
%Z = (X .* Y) + 2*Y.^2;
Z = log(sqrt(X .^ 2 + Y .^2));
end
function plotl(p1)
pi = 10*p1;
pj = -10*p1;
line([pi(1), pj(1)], [pi(2), pj(2)],'linewidth', 3, 'color','red');
end