clc;
clear;
close all;

M = [[314.63, 100];[100, 826.63]];
%M = eye(2);
y = [3.55; 2.24];
d = @(x) sqrt(dot(M*x,x));

%Barycentric Interpolation
cost = @(r, x) (d(x) - r*d(y))^2 + (d(y - x) - (1 - r)*d(y))^2;
cost = @(r, x)  r * d(x)/d(y) * d(x) + (1 - r)  * d(y - x)/d(y) * d(y - x);
figure
hold on
for r = 0.0:0.1:1.0
    h = @(x) cost(r, x);
    [res,c] = fminunc(h,[0;1]);
    %v = d(res) + d(y - res)
     c;
    plot(res(1), res(2),'o','color', 'blue','markersize',3);
    res2 = y*r;
    plot(res2(1), res2(2),'o','color', 'red');
    v = d(res2) + d(y - res2);
end
%projecao
%Na verdade, eu preciso adicionar na funcao de custo a propria restricao de
%x pertencer a reta. Isso pode ser feito pelo metodo da penalidade. 
%ref: https://gowerrobert.github.io/pdf/reports/GowerR_Painful_PCG_projections.pdf
proj_cost_real = @(x, p) d(p - x) + d(y - x);
proj_cost = @(x, p) d(p - x)^2 + d(y - x)^2; %adicionar restricao
p = [1;1];
h = @(x) proj_cost(x, p);
[res,c] = fminunc(h,[1;1]);
res;
proj = dot(p,y/norm(y)) * y/norm(y);

%gradient flow -> a way to obtain a geosedic when we only have a metric
figure 
hold on
cost = @(x, grad) d((x + grad/norm(grad)) - y)^2; % teria que adicionar uma restricao de d ser unitario ou fazer direto no S2
x = [0;0];
n = 25;
for i = 1:n
    h = @(grad) cost(x, grad);
    [res,c] = fminunc(h,[0;1]);
    x = x + (0.1)*res;
    plot(i/n * y(1), i/n * y(2),'o','color', 'red');
    plot(x(1), x(2),'o','color', 'blue','markersize',3);
end
keyboard