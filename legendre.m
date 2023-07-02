clc;
clear;
close all;

%legendre duality
p = [1;1];
x = linspace(-2,2,10);
l = [p(1),-p(2)];
y = p(1)*x - p(2);
y2 = 1*x + 2;
p2 = [1,-2];
figure
hold on
axis([-2,2,-2,2]);
%set(gca,'visible','off');
set(gcf,'color','white');
plot(p(1),p(2),'o','MarkerFaceColor','blue');
plot(x,y,'color','red');
plot(x,y2,'color','green');
plot(p2(1),p2(2),'o','MarkerFaceColor','cyan');
hold off

figure
hold on
axis([-2,2,-2,2]);
n = 6;
set(gcf,'color','white');
t = linspace(0,0+ 2*pi,n);
pts = [cos(t);sin(t)];
plot(pts(1,:),pts(2,:),'color','green');
plot(pts(1,:),pts(2,:),'o','MarkerFaceColor','blue');
hold off
figure
hold on
axis([-2,2,-2,2]);
for i = 1:n-1
    pt_i = pts(:,i);
    y_i = pt_i(1)*x - pt_i(2);
    pt_j = pts(:,mod(i,n-1) + 1);
    a = (pt_j(2) - pt_i(2)) / (pt_j(1) - pt_i(1));
    t = -pt_i(1)/(pt_j(1) - pt_i(1));
    p_corte = pt_i + (pt_j - pt_i)*t;
    b = -p_corte(2);
    y = pt_i(1)*x - pt_i(2);
    plot(a,b,'o','MarkerFaceColor','green');
    plot(x,y,'color','blue');

end

y = x.^2;
figure
hold on
axis([-2,2,-2,2]);
%set(gca,'visible','off');
set(gcf,'color','white');
plot(p(1),p(2),'o','MarkerFaceColor','blue');
plot(x,y,'color','red','o');
hold off