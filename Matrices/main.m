clear all;
close all;
clc;

P = [[1,0];[0,0]];
Q = [[0,1];[0,0]];
R = [[0,0];[1,0]];
S = [[0,0];[0,1]];
G = {P,Q,R,S};
for i = 1:4
    S*G{i};
end
d = pi/4;
D = [[2,0];[0,3]];
R = [[cos(d), -sin(d)];[sin(d), cos(d)]];
S = R' * D * R;
[V,D]  = eig(S);
D = diag(D);
dt = 0.000000001;
v = logm(3*eye(2));
l0 = v(1,1);
x = [1;1];
dt = 0.001;
sum = norm(x);
hs = [];
sum2 = 1;
for i  = 1:1/dt
    h = v*x*dt;
    x = x + h;
    sum = sum + norm(h);
    hs = [hs,h];
    sum2 = sum2 + (l0*dt)^i;
end
A = logm(S)
Q = V' * diag(log(D)) * V
T = (eye(2) + (A)*dt)^(1/dt)


