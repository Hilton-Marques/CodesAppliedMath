clc
clear all;
close all;

%Global variables
global dt KLR KLS u
% Grid Assemble
alpha = 1/100;
a = 0;
b = 1;
tf = 2;
n = 140;         % number of space gridpoints without boundaries
n2D = n^2;
dt = 0.1;                                          % time step                                        % final time
x = linspace(a,b,n+2)';
y = x;
xi = x(2:end-1);
h = x(2)-x(1);
lam = (alpha)*dt/(h^2);
disp(sprintf('CFL number: %0.2f',lam));
theta = 0.5;

%Matrix Assemble
I = speye(n);
I2D = speye(n2D);
e=ones(n,1);
A = spdiags([e,-2*e,e],-1:1,n,n);
A(1,1) = -1; %Neumann Boundary Condition
A(n,n) = -1; %Neumann Boundary Condition
A = kron(A,I) + kron(I,A);
KLR = (I2D - lam*theta*A); %Right Side
KLS = (I2D + lam*(1-theta)*A); %Left Side
% Initial Condition
u0 = f(xi,xi);
u = u0;

for tn = 1:ceil(tf/dt)
    %2ºQuestion - Linear problem
    %u = KLR\KLS*u;
    %u = gmres( KLR,KLS*u );
    %3ºQuestion - Non Linear problem
    xk = u;
    b = G(xk);
    u = JNFK(xk, b);
    clf
    % Plot 
    [X,Y] = meshgrid(xi);
    Z = griddata(reshape(X',[],1),reshape(Y',[],1),u,X,Y);
    Z0 = griddata(reshape(X',[],1),reshape(Y',[],1),u0,X,Y);
    surf(X,Y,Z);
    axis([0 1 0 1 0 1])
    title(sprintf('time t=%0.2f',tn*dt))
    drawnow
end

%Initial condition function
function z = f(x,y)
% initial condition function
[X,Y] = meshgrid(x,y);
Z = sin(pi*X).*sin(pi*Y);
surf(X,Y,Z)
z = reshape(Z',[],1);
end

% Source function
function s = S(x)
% Source function
s = 10*(x.^2 - x.^3);
end
% NonLinear function 
function G = G(xo)
global dt KLR KLS u
% Nonlinear function
G =  dt*0.5*S(xo) + dt*0.5*S(u)...
    - KLR*xo + KLS*u;
end
% Newtons Method
function x = JNFK(xk,b)
itermax = 50;        %Number of max iterations for Newton Mehthod         
tol = 10^-6;         % Tolerance for Newton Method
iter = 1;            % First iteration 

xn = xk - gmresJNFK(xk,b);  %First result
while (norm(xn - xk) > tol)
xk = xn;
b = G(xk);
xn = xk - gmresJNFK(xk,b);
iter = iter + 1;
if iter > itermax
    disp(sprintf('JNFK does not converge'));
    break;
end
end
x = xn;
end
% GMRES Jacobian Free
function x = gmresJNFK(xk,b)

restrt = 10;                    % Number of inner interations
max_it = 1;                     % Number of max interations
tol = 10^-6;                    % Gmres tolerance
ksi = 10^-8;                    % Tiny nudge for the Directional derivative
fxk = G(xk);                    % Function applied in the initial point

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
    r = (b -(G(xk + ksi*x)-fxk)/ksi);
    Q(:,1) = r / norm(r);
    s = norm( r )*e1;
    
    for i = 1:m                                   % construct orthonormal
        w = (G(xk + ksi*(Q(:,i))) ...
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
    r = (b - (G(xk + ksi*x/norm(x)) - ...
        G(xk))/ksi);                            % compute residual
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
%Pure GMRES
function x = gmres( A,b )

restrt = 10;
max_it = 1;
tol = 10^-6;

iter = 0;                                         % initialization
flag = 0;

[n,n] = size(A);                                  % initialize workspace
m = restrt;
Q(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(n,1);
e1(1) = 1.0;
x = b;
for iter = 1:max_it                              % begin iteration
    r = (b -A*x);
    Q(:,1) = r / norm( r );
    s = norm( r )*e1;
    
    for i = 1:m                                   % construct orthonormal
        w = (A*Q(:,i));                         % basis using Gram-Schmidt
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
    r = ( b-A*x );                              % compute residual
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




