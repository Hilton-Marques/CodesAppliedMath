clc;
clear all;
close all;

% Grid Assemble
alpha = 1;
a = 0;
b = pi;
ti = 0;
tf = 2;
n = 3;         % number of space gridpoints without boundaries
dt = 1;                                          % time step                                        % final time
x = linspace(a,b,n+2)';
t = linspace(ti,tf,(tf-ti)/dt);
xi = x(2:end-1); %internal points
h = x(2)-x(1)
lam = (alpha)*dt/(h^2);
disp(sprintf('CFL number: %0.2f',lam));
thetaMethod = struct('exp',0,'imp',1,'crN',0.5);
theta = thetaMethod.crN;
%Matrix Assemble
I = speye(n);
e=ones(n,1);
A = spdiags([e,-2*e,e],-1:1,n,n);
I = eye(n);
A = toeplitz ([-2 1 zeros(1,n-2)]);
A(1,1) = -2; %Dirichlet Boundary Condition
A(n,n) = -2; %Dirichlet Boundary Condition
KLR = (I - lam*theta*A); %Right Side
KLS = (I + lam*(1-theta)*A); %Left Side
%Boundary Condition
g = (lam*(1-theta) + lam*theta)*[0;zeros(n-2,1);0];
% Initial Condition
u0 = f(x);
u = u0(2:end-1);
for tn = 1:ceil(tf/dt)
    u = KLR\( (KLS*u)+g );
    clf
    %plot(x,u0,'b:',xi,u,'r.-')
    %title(sprintf('time t=%0.2f',tn*dt))
    %drawnow
    U(tn,:) = u;
end

%Plot surf

[X,Y] = meshgrid(xi,t);
surf(X,Y,U);
title(sprintf(strcat('Theta = ',num2str(theta))));
xlabel('comprimento');
ylabel('tempo');
zlabel('temperatura');

%Error
tf = 1.0;
x_controle = pi/2;
for i =1:2
    t=0;
    if (i==1)
    tl = 0.01;
    elseif (i==2)
    tl = 0.005;
    end
    l = (tl/0.4)^0.5;
    h(i) = l;
    dt(i) = tl; 
    nk = ceil(tf/dt(i));
    n = ceil((pi/h(i))) - 1;
    x = linspace(a,b,n+2)';
    x(2) - x(1)
    t = linspace(0,tf,nk);
    xi = x(2:end-1); %internal points
    [X,T] = meshgrid(xi,t);
    u0 = f(x);
    u = u0(2:end-1);
    lam = dt(i)/(h(i)^2);
    disp(sprintf('CFL number: %0.2f',lam));  
    %Matrix Assemble
    I = speye(n);
    e=ones(n,1);
    A = spdiags([e,-2*e,e],-1:1,n,n);
    A(1,1) = -2; %Neumann Boundary Condition
    A(n,n) = -2; %Neumann Boundary Condition
    KLR = (I - lam*theta*A); %Right Side
    KLS = (I + lam*(1-theta)*A); %Left Side
    for tn = 1:ceil(tf/dt(i))
    u = KLR\( (KLS*u));
    U_aste(tn,:) = u;
    end
    uexact = exp(-(T)).*sin(X);
    e = (uexact - U_aste);
    Err(:,i) = abs(max( e(end,:) ));
    if (i>1)
        p(i) = log ( Err(i-1)/Err(i) ) / log( dt(i-1) /dt(i) )
    end
    clear U_aste;
end

function y = f(x)
[n,~] = size(x);
%y = 20*ones(n,1); 
y =  sin(x);
end


