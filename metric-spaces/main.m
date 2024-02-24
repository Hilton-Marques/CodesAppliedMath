clc;
clear;
close all;

x = 1/3;
figure
hold on;
for i = 1:100
    plot(x,0, 'o', 'MarkerSize', 4, 'color', 'black');
    %text(x,0, num2str(i), 'VerticalAlignment','bottom');
    x = f(x);
end

M = [[144.45, -96.6]; ...
     [-96.6, 505.54]];
Minv = inv(M);
d = [1.73; -1.45];
v = [3.56; -0.21];
len = sqrt(dot(d, M*d));
len = sqrt(dot(v, M*v));
d = d/len
v = v/len
sqrt(dot(d, Minv*d))

function res = f(x)
res = x*x + 1/4;
end


