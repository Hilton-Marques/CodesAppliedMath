% Clear memory
clc;
clear;
%close all;
%Input Data
a = 0;
b = 1;
n = 100;
[xx,yy]=meshgrid(a:1/(n+1):b,a:1/(n+1):b);
f = ones(n,n);
f = -exp(xx(2:n+1,2:n+1).*yy(2:n+1,2:n+1)).*(xx(2:n+1,2:n+1).^2 + yy(2:n+1,2:n+1).^2);
%Boundary Condition Dirichlet
Bdr = zeros(n,n);
Bdr(1,:) = Bdr(1,:) + exp(xx(n+2,2:n+1) .* yy(n+2,1).*ones(1,n));
Bdr(n,:) = Bdr(n,:) + exp(xx(1,2:n+1) .* yy(1,1).*ones(1,n));
Bdr(:,1) = Bdr(:,1) + exp(yy(2:n+1,1) .* xx(1,1).*ones(n,1));
Bdr(:,n) = Bdr(:,n) + exp(yy(2:n+1,n+2) .* xx(1,n+2).*ones(n,1));
% F modified
f = (1/(n+1))^2*f + Bdr; 
% Fast Poisson Solver
L = 2 * ones(1, n) - 2 * cos((1: n)*pi/(n + 1));
LM = ones(n, 1) * L + L' * ones(1, n); %eigenvalues
AM = DSTTransform2D(f,n);
UM = DSTTransform2D(AM./LM,n);
surf(xx(2:n+1,2:n+1),yy(2:n+1,2:n+1),UM)
surf(xx(2:n+1,2:n+1),yy(2:n+1,2:n+1),exp(xx(2:n+1,2:n+1).*yy(2:n+1,2:n+1)))
max(max(UM))
max(max(exp(xx(2:n+1,2:n+1).*yy(2:n+1,2:n+1))))


%Outra forma de fazer
% zz=poisson(xx,yy);
% max(max(zz))
% surf(xx,yy,zz)

function B = DSTTransform2D(A,n)
v=[zeros(1,n); A; zeros(n+1,n)]; z = fft(v); % Fourier transform of size M
DSTRight = (2/(n+1))^0.5*-imag(z(2:n+1,1:n)) ;          % Sine transform of size N
v =[zeros(1,n); DSTRight'; zeros(n+1,n)]; z = fft(v);
DSTLeft =((2/(n+1))^0.5* -imag(z(2:n+1,1:n)))'    ;
B = DSTLeft;
end

function u=poisson(x,y)   % choose evaluation points for u(x,y)
% x = and y = choose 1 point or more
N=39; u=zeros(size(x));   % size=1,1 for evaluation at 1 point
if nargin==1              % Only x coordinates so 1D problem
    for k=1:2:N             % Add N terms in the 1D sine series
        u = u + 2^2/pi^3/k^3*sin(k*pi*x);
    end
    % xx=0:.01:1;yy=poisson(xx);plot(xx,yy) to plot u in 1D (2D is below)
    
elseif nargin==2          % x and y coordinates so 2D problem
    for i=1:2:N             % - u_xx - u_yy = 1 in unit square
        for j=1:2:N           % Add N^2 terms in the 2D sine series
            u = u + 2^4/pi^4/(i*j)/(i^2+j^2)*sin(i*pi*x).*sin(j*pi*y);
        end
    end                   % 3D would have (i*j*k)/(i^2+j^2+k^2)
end
end
