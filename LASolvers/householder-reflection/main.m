clear all;
clc;
%rng('default');
%invariants
%Case 2
b = [1;2;3];
a = 2;
V = [[1;0;0], [a;1;0]];
res = inv((V'*V))*V'*b;
easier = [b(1); b(2)] + [-a * b(2); 0];
assert(norm(easier - res) < 1e-8);
%Case 3
b = [1;2;3];
a = 2;
lam = 3.6;
V = [[1;0;0], [a; lam;0]];
res = inv((V'*V))*V'*b;
s = (1 - lam)/lam;
easier = [b(1); b(2)] + [-a * b(2) * (1/lam); s * b(2)];
assert(norm(easier - res) < 1e-8);
%Test
a = 2;
lam = 3;
A = [[1, a]; [0, lam]];
b = [3 ;2];
A\b;
s = (1 - lam);
easier = [b(1); b(2)] + (1/lam) * [-a * b(2); (1 - lam) * b(2)];
%Case 5
b_first = [5; 6];
v1 = [1;4];
v2 = [3;6];
V = [v1,v2];
proj_p = dot(v2,v1)/norm(v1);
a = proj_p/ norm(v1);
lam = -sqrt(dot(v2,v2) - proj_p^2) / norm(v1);
T = [[1, a]; [0, lam]];
theta = atan2(v1(2), v1(1));
Q = [[cos(theta), -sin(theta)]; [sin(theta), cos(theta)]];
T = 1/norm(v1) * Q';
A = T * V;
b = T * b_first;
s = (1 - lam);
easier = [b(1); b(2)] + (1/lam) * [-a * b(2); (1 - lam) * b(2)]
gab  = [v1,v2] \ b_first
keyboard;

%we consider m > n
A = rand(4,4);
[Q, R] = qr(A);
H = R' * R;
H2 = A' * A;
keyboard;
house(A);
function [Q,R] = house(A)
[m, n] = size(A);
U = zeros(m, n);

r1 = [-norm(A(:,1)); 0 ; 0];
w1 = A(:,1) - r1;
H1 = Reflection(w1);
u2 = H1 * A(:,2);
proj = [0; u2(2); u2(3)];
h2u2 = [0; norm(proj); 0];
w2 = proj - h2u2;
H2 = Reflection(w2);
r2 = H2 * u2;
r3 = H2 * H1 * A(:, 3);
R = [r1,r2,r3];
Q = (H2 * H1)';
[Q1, R1] = qr(A);
assert(Q1 == Q)
keyboard
end

function res = Reflection(v)
v = v/norm(v);
dim = size(v,1);
res = eye(dim) - 2 * v * v';
end