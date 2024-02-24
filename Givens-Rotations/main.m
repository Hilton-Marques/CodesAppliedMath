clc;
clear;
close all;

A = [[1;1;1],0.5 * [-1;-1;1]];

QRGivensRotations(A);


% figure
% hold on
% axis equal;
% axis off;
% light('Position', [-1.5,1,1], 'Style', 'local');
% lights = camlight;
% camproj('perspective');
% axis([-1,1,-1,1,-1,1])
% view(30,30);
% 
% 
% arrow3([0 0 0], A(:,1)');
% arrow3([0 0 0], A(:,2)');
% arrow3([0 0 0], [1,0,0], 'r');
% arrow3([0 0 0], [0,1,0], 'b');
% arrow3([0 0 0], [0,0,1], 'g');
% [Q,R] = qrgivens(A);
% 
% figure;
% hold on
% axis off
% quiver3(0,0,0, A(1,1), A(2,1), A(3,1));
% quiver3(0,0,0, A(1,2), A(2,2), A(3,2));
% view(30,30);
% 
% 
% function [Q,R] = qrgivens(A)
%   [m,n] = size(A);
%   Q = eye(m);
%   R = A;
% 
%   for j = 1:n
%     for i = m:-1:(j+1)
%       G = eye(m);
%       [c,s] = givensrotation( R(i-1,j),R(i,j) );
%       G([i-1, i],[i-1, i]) = [c -s; s c];
%       [V,D] = eig(G);
%       R = G'*R;
%       Q = Q*G;
%     end
%   end
% end
% 
% function [c,s] = givensrotation(a,b)
%   if b == 0
%     c = 1;
%     s = 0;
%   else
%     if abs(b) > abs(a)
%       r = a / b;
%       s = 1 / sqrt(1 + r^2);
%       c = s*r;
%     else
%       r = b / a;
%       c = 1 / sqrt(1 + r^2);
%       s = c*r;
%     end
%   end
% end

