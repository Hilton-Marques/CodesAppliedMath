%%POD 2D
%This code aims to implement the POD-DEIM
%method in a two-dimensional diffusion problem.
%The equation beign solved is: 
%u_t - (1/100)(u_xx+u_uyy) = 10(u^2 - u^3)

clear
clc
close all


%% Flag for DEIM method
flag_DEIM = true;

%% Mesh
dx = 0.00625;         %step for space
dt = 0.01;            %step for time
theta = 0.5;          %number for theta method
x = 0:dx:1;
y = x;
t = 0:dt:2;
n = length(x);
T = length(t);
sigma = 1/100;      %thermal conductivity
lam = (sigma)*dt/(dx^2);    %CFL number
disp(sprintf('CFL number: %0.2f',lam));

%% Matrix Assemble
e = ones(n,1);
A = spdiags([e -2*e e],-1:1,n,n);
A(1,1) = dx^2;           %Dirichlet Boundary Condition
A(1,2) = 0;
A(n,n) = dx^2;
A(n,n-1) = 0;
A_xx = A/dx/dx;
I = speye(n);
A = kron(A,I) + kron(I,A);
I2D = speye(n^2);
KLR = (I2D - lam*theta*A);          %Right Side
KLS = (I2D + lam*(1-theta)*A);      %Left Side
% Initial Condition
u0 = f0(x,y);
%Initializating matrixes
m = size(u0,1);
U_FD = zeros(m,T-1);      %Numerical Solution matrix
F = zeros(m,T-1);         %Source matrix
U_FD(:,1) = u0;
tic
%% Solve G_non_linear = 0 with NewtonJacobianFree
for i = 1:length(t) - 1
    xk = U_FD(:,i);
    F(:,i) = f(U_FD(:,i));         %Matrix F for DEIM
    b = dt*F(:,i) - KLR*U_FD(:,i) + KLS*U_FD(:,i);
    G_non_linear = @(xo) dt*0.5*(f(xo) + F(:,i))...
        - KLR*xo + KLS*U_FD(:,i);
    U_FD(:,i+1) = JNFK(xk, b,G_non_linear);
end
toc

%% Build POD basis
%Snapshots
L = 20;                               %Number of snapshots
snap = U_FD(:,1:ceil((T-1)/L):end);   %Equally spaced snapshots
%snap = U_FD(:,1:L);                  %First snapshots
L = size(snap,2);
singular_numbers = svd(snap);
%figure(3),semilogy(singular_numbers,'-o');
s = 1;
Sig = singular_numbers(s);
% Find the singular number greater than O(dx^2 + dt^2)
while(Sig >= (dx^2+dt^2))
    s = s + 1;
    Sig = singular_numbers(s);
    if s > L
        break;
    end
end
%% Build POD Matrixes
[U,S,V] = svd(snap,'econ');
Psi = U(:,1:s);
sol_pod(:,1) = Psi'*U_FD(:,1);
A_POD = Psi'*A*Psi;
I_POD = speye(s);
KLR_POD = (I_POD - lam*theta*A_POD);
KLS_POD = (I_POD + lam*(1-theta)*A_POD);

%% Build Deim Basis
[U_F,S_F,V_F] = svd(F(:,1:ceil((T-1)/L):end),'econ');
singular_numbers_F = diag(S_F);
p = 1;
Sig_F = singular_numbers_F(p);
while(Sig_F >= (dx^2+dt^2))
    p = p + 1;
    Sig_F = singular_numbers_F(p);
    if p > L
        break;
    end
end
Omega = U_F(:,1:p);
[~,Indices(1,1)] = max(abs(Omega(:,1)));
e1 = [zeros(Indices(1,1) - 1,1); 1; zeros(m -Indices(1,1),1)];
%Initialization of Permutation Matrix
P = zeros(m,p);
P(:,1) = e1;
for i = 2:p
    %c = (P(:,1:i-1)'*Omega(:,1:i-1))\(P(:,1:i-1)'*Omega(:,i));
    c = (Omega(Indices(1,1:i-1),1:i-1))\(Omega(Indices(1,1:i-1),i));
    r = Omega(:,i) - Omega(:,1:i-1)*c;
    [~,Indices(i)] = max(abs(r));
    P(:,i) = [zeros(Indices(1,i)-1,1); 1; zeros(m - Indices(1,i),1)];
end
%DEIM = Psi'*Omega*inv(P'*Omega);
DEIM = Psi'*Omega*inv(Omega(Indices(1,1:i),:));
dt_pod = dt;
t_pod = 0:dt_pod:2;
tic
%% Solve G_non_linear_POD = 0, with Newton Jacobian Free method
for i = 1:length(t_pod)-1
    xk = sol_pod(:,i);
    if (flag_DEIM)
        b = dt*DEIM*f(Psi(Indices,:)*xk) - KLR_POD*sol_pod(:,i) + KLS_POD*sol_pod(:,i);
        G_non_linear = @(xo) dt*0.5*DEIM*( f(Psi(Indices,:)*xo) + ...
            f(Psi(Indices,:)*sol_pod(:,i))) - KLR_POD*xo + KLS_POD*sol_pod(:,i);
    else
        b = dt*Psi'*f(Psi*xk) - KLR_POD*sol_pod(:,i) + KLS_POD*sol_pod(:,i);
        G_non_linear = @(xo) dt*0.5*Psi'*( f(Psi*xo) + ...
            f(Psi*sol_pod(:,i))) - KLR_POD*xo + KLS_POD*sol_pod(:,i);
    end
    sol_pod(:,i+1) = JNFK(xk, b,G_non_linear);
end
toc

%% Error
U_POD = Psi*sol_pod;
e = U_POD(:,end) - U_FD(:,end);
Err = vecnorm(e);

%% Plot final data
for i = 1:length(t) - 1
    clf
    % Plot
    j = i;
    if (i == 1)
        j = 0;
    end
    if (mod((j*dt),0.25) == 0)
        [X,Y] = meshgrid(x);
        Z_FD = griddata(reshape(X',[],1),reshape(Y',[],1),U_FD(:,i),X,Y);
        Z_POD = griddata(reshape(X',[],1),reshape(Y',[],1),U_POD(:,i),X,Y);
        subplot(3,1,1)
        surf(X,Y,Z_FD);
        title(sprintf('Solution FD-theta = 0.5, time t=%0.2f',i*dt))
        zlim([0 1]);
        set(gca,'FontSize',12)
        subplot(3,1,2)
        surf(X,Y,Z_POD);
        title(sprintf('Solution POD-DEIM-FD, time t=%0.2f',i*dt))
        set(gca,'FontSize',12)
        zlim([0 1]);
        subplot(3,1,3)
        surf(X,Y,Z_POD-Z_FD);
        title(sprintf('Error, time t=%0.2f',i*dt))
        set(gca,'FontSize',12)
        set(gcf,'position',[300,30,500,650])
        drawnow
        pause(0.1);
    end
    
end
function z = f0(x,y)
% initial condition function
[X,Y] = meshgrid(x,y);
Z = sin(pi*X).*sin(pi*Y);
z = Z(:);
end
function f = f(x)
% Source function
f = 10*(x.^2 - x.^3);
%f = zeros(size(x,1),1);   %Homogeneous problem
end
% Newtons Method
function x = JNFK(xk,b,G_non)
itermax = 100;        %Number of max iterations for Newton Mehthod
tol = 10^-5;         % Tolerance for Newton Method
iter = 1;            % First iteration

xn = xk - gmresJNFK(xk,b,G_non);  %First result
while (norm(xn - xk) > tol)
    xk = xn;
    b = G_non(xk);
    xn = xk - gmresJNFK(xk,b,G_non);
    iter = iter + 1;
    if iter > itermax
        disp(sprintf('JNFK does not converge'));
        break;
    end
end
x = xn;
end
function x = gmresJNFK(xk,b,G_non)

restrt = 10;                    % Number of inner interations
max_it = 1;                     % Number of max interations
tol = 10^-5;                    % Gmres tolerance
ksi = 10^-5;                    % Tiny nudge for the Directional derivative
fxk = G_non(xk);                    % Function applied in the initial point

iter = 0;                                         % initialization
flag = 0;

[n,~] = size(b);                                  % initialize workspace
m = restrt;
Q(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(n,1);
e1(1) = 1.0;
x = b;                                           % First Guess

for iter = 1:max_it                              % begin iteration
    r = (b -(G_non(xk + ksi*x)-fxk)/ksi);
    Q(:,1) = r / norm(r);
    s = norm( r )*e1;
    
    for i = 1:m                                   % construct orthonormal
        w = (G_non(xk + ksi*(Q(:,i))) ...
            - fxk)/ksi;
        % basis using Gram-Schmidt
        for k = 1:i
            H(k,i)= w'*Q(:,k);
            w = w - H(k,i)*Q(:,k);
        end
        H(i+1,i) = norm( w );
        Q(:,i+1) = w / H(i+1,i);
        
        for k = 1:i-1                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = cs(i)*s(i);                        % approximate residual norm
        s(i+1) = -sn(i)*s(i);
        s(i)   = temp;
        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
        error  = abs(s(i+1)) / norm(b);
        if ( error <= tol )                        % update approximation
            y = H(1:i,1:i) \ s(1:i);                 % and exit
            x = x + Q(:,1:i)*y;
            break;
        end
    end
    
    if ( error <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    x = x + Q(:,1:m)*y;                            % update approximation
    r = (b - (G_non(xk + ksi*x/norm(x)) - ...
        G_non(xk))/ksi);                            % compute residual
    s(i+1) = norm(r);
    error = s(i+1) / norm(b);                        % check convergence
    if ( error <= tol ), break, end
end

if ( error > tol ) flag = 1; end                % converged

    function [ c, s ] = rotmat( a, b )
        
        %
        % Compute the Givens rotation matrix parameters for a and b.
        %
        if ( b == 0.0 )
            c = 1.0;
            s = 0.0;
        elseif ( abs(b) > abs(a) )
            temp = a / b;
            s = 1.0 / sqrt( 1.0 + temp^2 );
            c = temp * s;
        else
            temp = b / a;
            c = 1.0 / sqrt( 1.0 + temp^2 );
            s = temp * c;
        end
    end

end