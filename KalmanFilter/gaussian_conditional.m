close all;
clc;
clear all;


%define two real symbolic variables
syms x y sigma_x sigma_y a b mu_x

h = a*x + b;

x_part = (-1/(2*sigma_x^2))*(x - mu_x)^2;
y_part = (-1/(2*sigma_y^2))*(y - h)^2;
p = expand(x_part + y_part);

[coef_x,xx] = coeffs(p, x);
[coef_y,yy] = coeffs(p, y);
[cross,x_y] = coeffs(p, [y,x]);

a  = simplify(coef_x(1));
latex(simplifyFraction(a))
b  = simplify(coef_y(1));
latex(simplifyFraction(b))
c = simplify(coef_x(2));
d = simplify(coef_y(2));
e = cross(end);
g = simplify(cross(2));
latex(simplifyFraction(g))


p = a * x^2 + b* y^2 + g * x * y + c*x + d*y + e; 
m = -2*[[a,g/2];[g/2, b]]
inv(m)
b = latex(simplifyFraction(e))
latex(simplifyFraction(a))
keyboard

