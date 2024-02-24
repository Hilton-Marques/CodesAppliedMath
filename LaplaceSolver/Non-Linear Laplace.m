
% Clear memory
clc;
clear;
close all;
%Input Data
global n K2D
a = 0;
b = 1;
n = 40; 
x = a:1/(n+1):b;
xi = x(2:end-1);
[xx,yy]=meshgrid(x);
f = ones(n,n);
f = -exp(xx(2:n+1,2:n+1).*yy(2:n+1,2:n+1)).*(xx(2:n+1,2:n+1).^2 + yy(2:n+1,2:n+1).^2);
%Boundary Condition Dirichlet
Bdr = zeros(n,n);
Bdr(1,:) = Bdr(1,:) + exp(xx(n+2,2:n+1) .* yy(n+2,1).*ones(1,n));
Bdr(n,:) = Bdr(n,:) + exp(xx(1,2:n+1) .* yy(1,1).*ones(1,n));
Bdr(:,1) = Bdr(:,1) + exp(yy(2:n+1,1) .* xx(1,1).*ones(n,1));
Bdr(:,n) = Bdr(:,n) + exp(yy(2:n+1,n+2) .* xx(1,n+2).*ones(n,1));
% F modified
f = (1/(n+1))^2*f ; 
% 
I = eye(n);
K = toeplitz([-2 1 zeros(1,n-2)]);
K2D = kron(K,I) + kron(I,K);
xk = ones(n^2,1);
b = G(xk);
u = JNFK(xk,b);
[X,Y] = meshgrid(xi);
Z = griddata(reshape(X',[],1),reshape(Y',[],1),u,X,Y);
surf(X,Y,Z)

function G = G(x)
global K2D
G = S(x) + K2D*x;
end
function s = S(x)
global n 
s = (1/(n+1))^2*(100*x + 1*(x./(ones(n^2,1) + x)));
end
function x = JNFK(xk,b)
itermax = 30;        %Number of max iterations for Newton Mehthod         
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

end
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

