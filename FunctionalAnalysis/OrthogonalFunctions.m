clear
close(findall(0,'Type','figure'));
clc;

x = linspace(-40,40);
y = linspace(-40,40);
[X,Y] = meshgrid(x,y);
Z = X.^2 - Y.^2;
Z_ortho = 2*X.*Y;
figure
hold on
axis([-40,40,-40,40]);
%surf(X,Y,Z)
%contour(X,Y,Z);
%contour(X,Y,Z_ortho);
plot(Z,Z_ortho,'b-x',Z',Z_ortho','b-x')

