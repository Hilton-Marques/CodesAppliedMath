clear all;
clc;
close all;

addpath(genpath("../../my_libs/"));

%Normal equation
Cinv = [[8.5, 1.04;
			1.04, 4.64]]/8;
C = inv(Cinv)
H = [1;1];
b = [4;3];
H*(H'*Cinv*H)^-1*H'*Cinv*b

% Interpolation by Sylvester Law of Inertia
angle = pi/2;
D = [[cos(angle), -sin(angle)];[sin(angle), cos(angle)]];
E = [[1,0];[0,100]];
Q1 = D' * E * D;
angle = pi/3;
D = [[cos(angle), -sin(angle)];[sin(angle), cos(angle)]];
E = [[1,0];[0,200]];
Q2 = D' * E * D;

% Q1 = [[4,0];[0,4]];
% Q2 = [[51,-50];[-50,51]];

[V1,D1] = eig(Q1);
[V2,D2] = eig(Q2);
C = chol(Q2)'*inv(chol(Q1)');
Q2_new = C*Q1*C'
check = norm(Q2_new - Q2);

q = QuadraticForms(Q1);
A1 = q.Translate(Q2,[1;0]);
A3 = q.TranslateEuclidean(Q2,[1;2]);
A2 = q.TranslateRiemmannian(Q2,[1;4]);
q.exportFrame(filename="euclidean_and_cholesky_rieman");

