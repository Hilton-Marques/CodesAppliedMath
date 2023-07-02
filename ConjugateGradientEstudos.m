clc;
clear all;
close all;
%Para entender o conjugate gradient
%Quadratic form
A = [5 4;4 5];
b = [2 -8]';
x = -1000:100:1000;
[X,Y] = meshgrid(x);
P(1,:) = reshape(X,1,[]);
P(2,:) = reshape(Y,1,[]);
f = 0.5*P'*A*P - b'*P;
z = diag(f)';
Z = griddata(P(1,:),P(2,:),z,X,Y);
surfc(X,Y,Z);
[U,V] = gradient(Z,1,1);
hold on
quiver(X,Y,U,V)
hold off
%Como encontrar a ellipse
A= [5 4 ;4 5] % You can change A
V = 1/sqrt(2)*[1 1;1 -1]; %eigenvectors
D = [9 0;0 1];
nD = [1/3 0;0 1];%Os eixos da ellipse
H = V*nD*V'; %Essa matriz transforma um circulo na ellipse que eu quero
theta= [0:2*pi/50 :2*pi]; %circulo
circle= [cos(theta); sin(theta)];
ellipse = H * circle; %ellipse que eu quero
diag(ellipse'*A*ellipse); %Os pontos na ellipse tem essa propriedade
axis([-4 4 -4 4]); axis('square')
%As duas ellipses com essa propriedade, muito lindo!
%O grande challenge, transformar a ellipse em um círculo
L = chol(A);
B = (D.^0.5)*V';
circle = B*ellipse;
circle = L*ellipse;
figure,plot(circle(1,:), circle(2,:),ellipse(1,:),ellipse(2,:))
%Outra forma de obter a ellipse
theta= [0:2*pi/50 :2*pi]; %circulo
circle= [cos(theta); sin(theta)];
nellipse = inv(L)*circle;
figure,plot(circle(1,:), circle(2,:),nellipse(1,:),nellipse(2,:));

  