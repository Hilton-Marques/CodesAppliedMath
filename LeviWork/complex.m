clc;
clear;
close all;

figure
hold on
axis equal
axis off
axis([-1,5,-2,2]);
n = 30;
teta = linspace(0,2*pi,n);
z = exp(complex(0,teta));
plot(real(z),imag(z));
plotLine(0,6);
proj = ray2([0.2,2]);
%plot(real(proj),imag(proj),'color','black');
c = orthoCircle(0,1);
%plot(real(c),imag(c));
homotopy(linspace(0,10,100))
h1 = drawCurve(orthoCircle(1/2,1));
%rays(ray3(orthoCircle(1/5,1)));
%ray(orthoCircle(1/5,1))
h2 = drawCurve(orthoCircle(1,2));
len = lenCurve(ray3(orthoCircle(1,2)));
homotopy(orthoCircle(1/5,1))
delete(h1);
homotopy(orthoCircle(1,5))
delete(h2);
function plotLine(p1,p2)
line([real(p1),real(p2)],[imag(p1),imag(p2)],'color','black');
end
function proj = ray(pt)
N = complex(0,1);
u = pt - N;
proj = N + 2*u;
for i = 1:length(pt)
    line([real(N),real(proj(i))],[imag(N),imag(proj(i))],'color','yellow');
end
end
function len = lenCurve(curve)
len = 0;
for i = 1:length(curve) - 1
    pi = curve(i);
    pj = curve(i + 1);
    len = len + norm(pj - pi);
end
end
function proj = ray2(pt)
N = complex(0,1);
R = sqrt(2);
u = pt - N;
proj = N + (R*R)./conj(u);
end
function homotopy(curve)
h = plot(real(curve),imag(curve),'color','black');
pause(1);

final = ray3(curve);
n = length(curve);
delete(h);
t = 50;
for i = 1:t
    inter = curve + (i/t)*(final - curve);
    h = plot(real(inter),imag(inter),'color','black');
    pause(0.01);
    delete(h);
end
plot(real(final),imag(final),'color','black');
end
function proj = ray3(pt)
i = complex(0,1);
proj = (1 + i *pt) ./ (i + pt);
end
function c = orthoCircle(p0,p1)
avg = (p0 + p1)/2;
R = abs(avg - p1);
theta = linspace(0,pi,1000000);
c = R*exp(complex(0,theta)) + avg;
end
function h = drawCurve(curve)
h = plot(real(curve),imag(curve),'color','black');
end
function rays(pt)
N = complex(0,1);
for i = 1:length(pt)
    u = pt(i) - N;
    scale = norm(u);
    v = 2*u/scale + N;
    line([real(N),scale*real(v)],[imag(N),scale*imag(v)],'color','yellow');
end
end