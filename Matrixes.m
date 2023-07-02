clc
n = 2;
k=2;
t = pi;
[K,T,B,C] = KTBC(n);
x = (1:n)'.^2; % squares 
x = [zeros(k-1,1) ; 1 ; zeros(n-k,1)]; % delta
%x = [zeros(k-1,1) ; ones(n-k+1,1)]; % step
%x = [zeros(k-1,1) ; (0:(n-k))']; % ramp
%x = exp((1:n)'*i*t)
mesh = linspace(0,1,n+2)
x = sin(mesh(2:n+1)*t)
[V,D]=eig(1000*K);
D
A=chol(K)
KKT = [1 0 -1; 0 1 1; -1 1 0];
[L,U]= lu(KKT);
L'
diag(U)
function [K,T,B,C] = KTBC(n)
    % Create the four special matrices assuming n>1
    K = toeplitz ([2 -1 zeros(1,n-2)]);
    T = K; T(1,1) = 1;
    B = K; B(1,1) = 1; B(n,n) = 1;
    C = K; C(1,n) = -1; C(n,1) = -1;
end
