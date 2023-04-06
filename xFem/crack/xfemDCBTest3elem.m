close all;
clear;
clc;

set(gcf,'color','w');

d1 = 10;
d2 = 10;
d3 = -10;
d4 = -10;

nodes = [ 0.00,  0.00;
20.00,  0.00;
40.00,  0.00;
60.00,  0.00;
 0.00,  20.00;
20.00,  20.00;
40.00,  20.00;
60.00,  20.00];

conec = [1,  2,  6,  5;
2,  3,  7,  6;
3,  4,  8,  7];

interface = [61, 10.0,0, 10.0];

u = [   0.0000,  0.0000
       -0.0000, -0.0125;
       -0.0000, -0.0250;
       -0.0000, -0.0375;
        0.0000,  0.0000;
        0.0000,  0.0125;
        0.0000,  0.0250;
        0.0000,  0.0375];

a = [-0.0063, -0.0000;
     -0.0063, -0.0125;
     -0.0063, -0.0250;
     -0.0063, -0.0375;
      0.0062,  0.0000;
      0.0062, -0.0125;
      0.0062, -0.0250;
      0.0063, -0.0375];

n = 50;
qsi = zeros(n, n);
eta = zeros(n, n);
for i = 1:n
    qsi (:,i) = 2/(n - 1)*i + (-1 - 2/(n - 1));
    eta (i,:) = - 2/(n - 1)*i + (1 + 2/(n - 1));
end
[TX,TY] = meshgrid(1:n);
%Connectivity Matrix
c = 0;
conect = zeros((size(TX,1)-1)*(size(TX,2)-1),4);
for j = 0:(size(TX,2)-2)
    for i = 1: size(TX,1)-1
        c = c + 1;
        ii = j*(size(TX,1)) + i;
        conect(c,:) = [ii, ii + 1, ii + size(TX,1) + 1, ii + size(TX,1)];
    end
end
N1std = 0.25*(1-qsi).*(1-eta);
N2std = 0.25*(1+qsi).*(1-eta);
N3std = 0.25*(1+qsi).*(1+eta);
N4std = 0.25*(1-qsi).*(1+eta);

heav1 = heaviside(n, qsi, eta, d1, d2, d3, d4);

%Campo de enriquecimeto 1
N1enr = N1std.*(heav1 - heav1(n,1)*ones(n, n));
N2enr = N2std.*(heav1 - heav1(n,n)*ones(n, n));
N3enr = N3std.*(heav1 - heav1(1,n)*ones(n, n));
N4enr = N4std.*(heav1 - heav1(1,1)*ones(n, n));

%Campo de deslocamento
xini = qsi*10;
yini = eta*10;
%figure (1)
hold on
axis off
%grid on
for i = 1:size(interface,1)
    crack = interface(i,:);
    %line([crack(1),crack(3)],[crack(2),crack(4)],'linewidth',2,'color','cyan')
end
squares = Square();
for i = 1:size(conec,1)
 %coordenadas iniciais
  xe = (N1std*nodes(conec(i,1),1) + N2std*nodes(conec(i,2),1) + N3std*nodes(conec(i,3),1) + N4std*nodes(conec(i,4),1));
  ye = (N1std*nodes(conec(i,1),2) + N2std*nodes(conec(i,2),2) + N3std*nodes(conec(i,3),2) + N4std*nodes(conec(i,4),2));
 %campo de deslocamento
  uxe = (N1std*u(conec(i,1),1) + N1enr*a(conec(i,1),1) + N2std*u(conec(i,2),1) + N2enr*a(conec(i,2),1) + N3std*u(conec(i,3),1) + N3enr*a(conec(i,3),1) + N4std*u(conec(i,4),1) + N4enr*a(conec(i,4),1));
  uye = (N1std*u(conec(i,1),2) + N1enr*a(conec(i,1),2) + N2std*u(conec(i,2),2) + N2enr*a(conec(i,2),2) + N3std*u(conec(i,3),2) + N3enr*a(conec(i,3),2) + N4std*u(conec(i,4),2) + N4enr*a(conec(i,4),2));
 %elemento deformado
 xed = xe + 100*uxe;
 yed = ye + 100*uye;
 s = sqrt(uxe .^ 2 + uye.^2);
 pcolor(xed,yed,s);
 shading interp;
 colormap('jet');
%plot(xed,yed,'o','color','black');
 s = s(:);

%DT = delaunay(xed(:),yed(:));
%triplot(DT,xed(:),yed(:));
pts = [xed(:),yed(:)];
for j = 1:size(conect,1)
    sq = Square(pts(conect(j,:),:),s(conect(j,:)));
    if sq.crackIntersection(interface)
        sq.show();
        squares(end+1) = sq;
    end
end
 %plot(xed,yed,'o','color','b');
end
%plot mesh
figure (2)
hold on
grid on
for i = 1:size(conec,1)
  xe = [nodes(conec(i,1),1); nodes(conec(i,2),1); nodes(conec(i,3),1); nodes(conec(i,4),1); nodes(conec(i,1),1)];
  ye = [nodes(conec(i,1),2); nodes(conec(i,2),2); nodes(conec(i,3),2); nodes(conec(i,4),2); nodes(conec(i,1),2)];
  plot(xe, ye)
end
for i = 1:(size(interface,1))
  ix = [interface(i,1); interface(i,3)];
  iy = [interface(i,2); interface(i,4)];
  plot(ix, iy, "-b")
end
set(gca,'FontSize',15)

%plot displacement + mesh
figure (3)
hold on
grid on
for i = 1:size(conec,1)
  xe = [nodes(conec(i,1),1); nodes(conec(i,2),1); nodes(conec(i,3),1); nodes(conec(i,4),1); nodes(conec(i,1),1)];
  ye = [nodes(conec(i,1),2); nodes(conec(i,2),2); nodes(conec(i,3),2); nodes(conec(i,4),2); nodes(conec(i,1),2)];
  plot(xe, ye)
  xed = [nodes(conec(i,1),1) + u(conec(i,1),1); nodes(conec(i,2),1) + u(conec(i,2),1); nodes(conec(i,3),1) + u(conec(i,3),1); nodes(conec(i,4),1) + u(conec(i,4),1); nodes(conec(i,1),1) + u(conec(i,1),1)];
  yed = [nodes(conec(i,1),2) + u(conec(i,1),2); nodes(conec(i,2),2) + u(conec(i,2),2); nodes(conec(i,3),2) + u(conec(i,3),2); nodes(conec(i,4),2) + u(conec(i,4),2); nodes(conec(i,1),2) + u(conec(i,1),2)];
  plot(xed, yed, "-r")
end
for i = 1:(size(interface,1))
  ix = [interface(i,1); interface(i,3)];
  iy = [interface(i,2); interface(i,4)];
  plot(ix, iy, "-b")
end
set(gca,'FontSize',15)
