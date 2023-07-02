% Clear memory
clc;
clear;
%close all;

%Input data

a = 0;  % begin of interval
b = 6;  % end of interval
c = 1;  % constant
Ti = 0; %Dirichlet
Tf = 5; %Dirichlet
Neu = 1;
dT = -5; %Neumann
h = 1.5;
n = ((b-a)/h -1); % size of interval
mesh = linspace(a,b,n+2)';
if (Neu == 1)
   n = n+1;
end 
f = (h^2)*-0.247*(3*ones(n,1) - mesh(2:n+1,1)).^2; % f function on internal nodes
% Solver
K = toeplitz ([-2 1 zeros(1,n-2)]);
%Boundary Condition
K(1,1) = -2;
K(n,n) = (-2 + 2*h);
K(n,n-1) = 2;
f(1,1) = f(1,1) - (Ti);
f(n,1) = f(n,1) - 2*h*(dT);
u = K\f;

