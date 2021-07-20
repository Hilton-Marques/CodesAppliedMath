clc;
clear all;
close all;
% Subspace Kn
n = 4; %Subspace de Kyrlov
A = diag([1 2 3 4]);
A = rand(n);
xo = [1 1 1 1]';
%Arnlod's Orthogonalization
H = zeros(4,n);
Q = zeros(4,n);
q1 = xo/vecnorm(xo);
Q(:,1) = q1;
for j =1:n-1
    t = A*Q(:,j);
    for i=1:j
        H(i,j)=Q(:,i)'*t;
        t = t - H(i,j)*Q(:,i); %Garante que todos os q's sao ortogonais
    end
    H(j+1,j) = vecnorm(t);
    Q(:,j+1) = t/H(j+1,j);
end
H(:,n) = Q'*A*Q(:,n)
Q;
H;

