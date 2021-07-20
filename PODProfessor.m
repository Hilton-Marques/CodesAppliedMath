clear
clc
close all

dx = 0.0125;
%CFL = 0.3;
dt_snap =dx;
dt_pod = dx;

T = 1;
x = 0:dx:1;
t=  0:dt_snap:T;
n = length(x);
sigma = 1/10;


e = ones(n,1);

 A = spdiags([e -2*e e],-1:1,n,n);
 A(1,1) = dx^2;
 A(1,2) = 0;
 A(n,n) = dx^2;
 A(n,n-1) = 0;
 A_xx = A/dx/dx;

count = 0;
figure
%for sigma = [0 1/2 1/10 1/100]
    count = count +1;
    A = -A_xx;
    I = speye(n);
    u(:,1,count) = max(sin(pi*x),0);

tic
for i = 1:length(t)-1
    u(:,i+1,count) = (I+dt_snap*sigma*A)\u(:,i,count);
end
%valores{count} = svd(u)
    tmp = svd(u(:,:,count));
     semilogy(tmp,'*','LineWidth',1)
     
     hold on
toc
%end
grid
%legend('0' ,'1/2', '1/10', '1/100')
set(gca,'FontSize',16)


figure
surf(x,t,u')
axis([0 x(end) 0 T -0.2 1])
xlabel('x')
ylabel('t')
title('Solution of heat equation')
set(gca,'FontSize',16)

%% Proper Orthogonal Decomposition

% snapshosts
snap = u(:,1:4:end); %matriz dos snapshots

[U,S,V] = svd(snap);
valores = diag(S);
figure
semilogy(svd(snap),'*-','LineWidth',2)
grid
title('Decay Singular Values')
%keyboard
ratio = 0;
l= 0
% achamos tamanho problema reduzido
while(ratio<0.99)
l= l+1;
ratio = (sum(valores(1:l).^2)/sum(valores.^2));
end

%semilogy(1:l,ratio,'LineWidth',2)
%grid



%for l = 1:rank(snap)

Psi = U(:,1:l); % base POD

sol_pod(:,1) = Psi'*u(:,1); %projeção do dado inicial
A_pod = Psi'*A*Psi; %projeção da matriz A
I_pod = speye(l);

t_pod = 0:dt_pod:T;
tic
for i =1:length(t_pod)-1
    sol_pod(:,i+1) = (I_pod+dt_pod*sigma*A_pod)\sol_pod(:,i); 
end
time = toc
%error(l) = max(max(abs(u' -(Psi*sol_pod)')));
%clear sol_pod
%end
%%
% figure
% semilogy(1:l,error,'LineWidth',3)
% title('Inf error POD')
% grid
% set(gca,'FontSize',16);

%%

% 
figure
surf(x,t_pod,(Psi*sol_pod)')
axis([0 x(end) 0 T -0.2 1.1])
xlabel('x')
ylabel('t')
title('Solution of POD heat equation')
set(gca,'FontSize',16)
% 
figure
surf(x,t_pod,abs(u' -(Psi*sol_pod)'))
xlim([0 x(end)])
ylim([0 T])
%axis([0 1 0 2 0 1])
xlabel('x')
ylabel('t')
title('Absolute difference betwenn Solution and POD')
set(gca,'FontSize',16)


