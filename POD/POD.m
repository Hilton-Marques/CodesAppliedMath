clear
clc
close all


dx = 0.01;
dt = 0.025;
L = 5; % snapshots

dt_snap = 0.025;
dt_pod = 0.025;

x = 0:dx:1;
t = 0:dt_snap:2;

n = length(x);
sigma = 1/10;
e = ones(n,1);
A = spdiags([e -2*e e],-1:1,n,n);
A(1,1) = dx^2;
A(1,2) = 0;
A(n,n) = dx^2;
A(n,n-1) = 0;
A_xx = A/dx/dx;
I = speye(n);
u(:,1) = max(sin(pi*x),0);
tic
for i = 1:length(t) - 1
    u(:,i+1) = (I - dt*sigma*A)\u(:,1);
end
figure(1),surf(x,t,u');
toc
%% Proper Othogonal

%snapshots
snap = u(:,1:L); %matriz dos snapshots
figure(3),semilogy(svd(snap),'-o');
[Un,Sn] = eig(snap'*snap);
[U,S,V] = svd(snap);
valores = diag(S);
ratio = 0;
M = 0;
% achamos tamanho problema reduzido
while(ratio<0.99)
    M = M +1;
    ratio = (sum(valores(1:M).^2)/sum(valores.^2));
end

Psi = U(:,1:M);
sol_pod(:,1:L) = Psi'*u(:,1:L);

A_pod = Psi'*A*Psi;
I_pod = speye(M);
%t_pod = 0:dt_pod:2;
t_pod = dt_snap*L:dt_snap:2;
tic
for i = 1:length(t_pod)
    sol_pod(:,L+i) = (I_pod - (dt_snap)*sigma*A_pod)\ ...
        sol_pod(:,L-1+i);
end
toc
Uaprox = Psi*sol_pod;
e = Uaprox - u;
e = vecnorm(e);

figure(2),surf(x,t,(Psi*sol_pod)')
