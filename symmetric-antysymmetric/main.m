%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The symmetric and antysemetric parts of matrix induces a  
% decomposition of the space. However this decomposition is not orthogonal
%as see in the dot product below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

A = [[1,2];[3,4]];

Ap = (A + eye(2))/2;
Aq = (A - eye(2))/2;
Ap + Aq
x = [1;1];
y = A*x

figure;
hold on;
y1 = Ap * x;
y2 = Aq * x;
y1 + y2
quiver(0, 0, y1(1), y1(2));
quiver(0, 0, y2(1), y2(2));
quiver(0, 0, y(1), y(2));
dot(y1,y2)
