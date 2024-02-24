clear all;
close all;
clc;

syms x y lam;

f(x,y) = x^7 + y^7;
d = [-0.5000; -0.3333];   
x0 = [3;2];
g = x0 + lam*d;
h = f(g(1),g(2));
g = diff(h, lam)
solve(g == 0, lam, 'MaxDegree', 4);
extrema = vpa(ans, 6)

fplot(h)
hold on
plot(extrema, subs(h,extrema), '*')
hold off
