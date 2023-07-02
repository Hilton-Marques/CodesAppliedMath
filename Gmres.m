clc;
clear all;
close all;
A = randi(10, [4,4]);

A = [8     7     6     1;...
     6     7     2     5;...
     5     9     3     2;...
    10     9     9    10];
b = randi([0, 1], [4,1]);
 b= [1;0;1;0];
x = gmres(A,b)
A\b
function x = gmress( A,b )

restrt = 10;
max_it = 1;
tol = 10^-6;
M = diag(A);
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
            y = H(1:i,1:i)\ s(1:i);                 % and exit
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