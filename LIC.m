clc;
clear all;
close all;
set(gcf,'color','w');
set(gca,'visible','off');
set(gca,'xtick',[]);
step = 0.1;
x = 0:step:1;
[X,Y]=meshgrid(x);
P = [X Y];
Q = [Y X];
hold on
plot(P,Q,'Color','black');

for j = 1:length(x)-1
    for i= 1:length(x)-1
        xplot = x(1,i) + step/2;
        yplot = x(1,j) + step/2;
        xi = 2*xplot - 1;
        yi = 2*yplot - 1;
        Vx = -yi;
        Vy = xi;
        v = [Vx;Vy];
        v = (0.07)*v/vecnorm(v);
        plot(xplot,yplot,'o','Color','red');
        xc = [xplot;yplot];
        vector2D(v+xc,[0 0 1],xc);
        axis([0 1 0 1]);
    end
end
for h=1:2
k = 7;
l = 4+h;
StepLen = 0.1;
L = 10;
xplot = x(1,k) + step/2;
yplot = x(1,l) + step/2;
xc = [xplot;yplot];
lineplot = zeros(2,L);
lineplot(:,1)=xc;
%Highlight fist pixel
c = [0 1 0];
if (h==2)
 c = [1 0 0];
end
fill([x(1,k) x(1,k+1) x(1,k+1) x(1,k)],...
    [x(1,l) x(1,l) x(1,l+1) x(1,l+1)],...
    c,'facealpha',0.5);


for i= 2:L
    xplot = x(1,k) + step/2;
    yplot = x(1,l) + step/2;
    xi = 2*xplot - 1;
    yi = 2*yplot - 1;
    Vx = -yi;
    Vy = xi;
    v = [Vx;Vy];
    v = v/vecnorm(v);
    lineplot(:,i) = lineplot(:,i-1) + StepLen*v;
    linef = lineplot(:,i);
    linei = lineplot(:,i-1);
    c = [0 1 0];
    if (h==2)
        c = [1 0 0];
    end
    line([linei(1,1) linef(1,1)],[linei(2,1) linef(2,1)],'LineWidth',2,'Color',c);
    axis([0 1 0 1]);
    drawnow;
    k = ceil(linef(1,1)/step);
    l = ceil(linef(2,1)/step);
    pause(1);
end
end
keyboard
hold off;
clf;
%%LIC initial Data
k = 7;
l = 6;
StepLen = 0.1;
L = 10;
xplot = x(1,k) + step/2;
yplot = x(1,l) + step/2;
xc = [xplot;yplot];
lineplot = zeros(2,L);
lineplot(:,1)=xc;
noise = zeros(10,10,3);
for m = 1:L
    for n= 1:L
        noise(m,n,:) = rand*[1 1 1];
    end
end
%%plot
set(gcf,'color','w');
subplot(2,4,8);
image(noise);
set(gca,'visible','off');
set(gca,'xtick',[]);
subplot(2,4,[1,2,3,5,6,7]);
set(gca,'visible','off');
set(gca,'xtick',[]);
hold on;
plot(P,Q,'Color','black');
%Highlight fist pixel
fill([x(1,k)*(1-0.01) x(1,k)*1.01 x(1,k)*1.01 x(1,k)*(1-0.01)],...
     [x(1,l) x(1,l) x(1,l+1) x(1,l+1)], ...
     [0.5 0.5 0],...
     [x(1,k+1)*(1-0.01) x(1,k+1)*1.01 x(1,k+1)*1.01 x(1,k+1)*(1-0.01)],...
     [x(1,l) x(1,l) x(1,l+1) x(1,l+1)], ...
     [0.5 0.5 0],...
     [x(1,k) x(1,k+1) x(1,k+1) x(1,k)],...
     [x(1,l)*(1-0.01) x(1,l)*(1-0.01) x(1,l)*1.01 x(1,l)*1.01], ...
     [0.5 0.5 0],...
     [x(1,k) x(1,k+1) x(1,k+1) x(1,k)],...
     [x(1,l+1)*(1-0.01) x(1,l+1)*(1-0.01) x(1,l+1)*1.01 x(1,l+1)*1.01],...
     [0.5 0.5 0]);
subplot(2,4,4);
set(gca,'visible','off');
set(gca,'xtick',[]);
hold on;
x_F = linspace(-L,L,50);
F = 0.25*(1 + cos((pi/L)*x_F)).*(1+cos(x_F + 0));
plot(x_F,F);
for i= 2:L
    subplot(2,4,[1,2,3,5,6,7]);
    xplot = x(1,k) + step/2;
    yplot = x(1,l) + step/2;
    xi = 2*xplot - 1;
    yi = 2*yplot - 1;
    Vx = -yi;
    Vy = xi;
    v = [Vx;Vy];
    v = v/vecnorm(v);
    xc = [xplot;yplot];
    lineplot(:,i) = lineplot(:,i-1) + StepLen*v;
    h1 = fill([x(1,k) x(1,k+1) x(1,k+1) x(1,k)],[x(1,l) x(1,l) x(1,l+1) x(1,l+1)],...
        0.5*[i/L i/L i/L],'facealpha',0.5);
    plot(xplot,yplot,'o','Color','red');
    %vector2D(0.05*v+xc,[0 0 1],xc);
    linef = lineplot(:,i);
    linei = lineplot(:,i-1);
    %line([linei(1,1) linef(1,1)],[linei(2,1) linef(2,1)],'LineWidth',2,'Color','gree');
    axis([0 1 0 1]);
    drawnow;
    k = ceil(linef(1,1)/step);
    l = ceil(linef(2,1)/step);
    subplot(2,4,4);
    h2 = fill( [i-2 (i-2)*1.01 (i-2)*1.01 i-2] , ...
    [0 0 filter((i-2),L) filter(((i-2)*1.01),L)],[ 0 1 0]);
    axis([-L L 0 1]);
    pause(1);
    delete(h2);

end
L = 10;
k = 7;
l = 6;
xplot = x(1,k) + step/2;
yplot = x(1,l) + step/2;
xc = [xplot;yplot];
lineplot = zeros(2,L);
lineplot(:,1) = xc;
for i= 2:L
    subplot(2,4,[1,2,3,5,6,7]);
    xplot = x(1,k) + step/2;
    yplot = x(1,l) + step/2;
    xi = 2*xplot - 1;
    yi = 2*yplot - 1;
    Vx = -yi;
    Vy = xi;
    v = [Vx;Vy];
    v = -v/vecnorm(v);
    xc = [xplot;yplot];
    lineplot(:,i) = lineplot(:,i-1) + StepLen*v;
    h1 = fill([x(1,k) x(1,k+1) x(1,k+1) x(1,k)],[x(1,l) x(1,l) x(1,l+1) x(1,l+1)],...
        0.5*[i/L i/L i/L],'facealpha',0.5);
    plot(xplot,yplot,'o','Color','red');
    vector2D(0.05*v+xc,[0 0 1],xc);
    linef = lineplot(:,i);
    linei = lineplot(:,i-1);
    line([linei(1,1) linef(1,1)],[linei(2,1) linef(2,1)],'LineWidth',2,'Color','red');
    axis([0 1 0 1]);
    drawnow;
    k = ceil(linef(1,1)/step);
    l = ceil(linef(2,1)/step);
    subplot(2,4,4);
    h2 = fill( [-(i-2) -(i-2)*1.01 -(i-2)*1.01 -(i-2)] , ...
    [0 0 filter((i-2),L) filter(((i-2)*1.01),L)],[ 0 1 0]);
    axis([-L L 0 1]);
    pause(1);
    delete(h2);
end
function triangle(xf,c,xi)
x = xf - xi;
h = vecnorm(x)*0.25;
b = vecnorm(x)*0.25;

p_x = [-x(2,1); x(1,1)];
p_x = p_x/norm(p_x);
p1 = xf + b/2*p_x;
p2 = xf - b/2*p_x;
p3 = xf + h*(x/norm(x));

X = [p1(1,1), p2(1,1) , p3(1,1)];
Y = [p1(2,1), p2(2,1) , p3(2,1)];
Thandle = fill(X, Y, c);
end

function vector2D(xf,c,xi)
if (nargin < 2)
    c = [1 0 1];
end
if (nargin < 3)
    xi = zeros(2,1);
end
X = [xi(1,1), xf(1,1)];
Y = [xi(2,1), xf(2,1)];
line( X, Y, 'Color', c,'Linewidth',2);
triangle(xf,c,xi);
end
function f = filter(x,L)
f = (0.25*(1 + cos((pi/L)*x))*(1+cos(x + 0)));
end