close all;
clc;
clear all;

addpath(genpath("../../my_libs/"));
addpath(genpath("../../RVC3-MATLAB/toolbox"));

a = SO3().getRandomElement
b = Manifold.skew([3,2,1]);
c = a*b;
d = b*a;
L = 25/100;     % m
w = 18/100;     % m
b = 2.5/100;    % m
m = 0.7;        % kg
w0 = [0.1, 20, 0.1]';
dt = 0.005;
obj = RotationalMechanics(L,w,b,m,w0,dt);
obj.Solver();

scene = Scene('geodesics',"SE3");
% Cube element

g = 0;        % m/s^2
scale = diag([L,w,b]);
% Calculate the principal moments of inertia:

lambda1 = (m/12)*(w^2 + b^2); 	% kg-m^2
lambda2 = (m/12)*(L^2 + b^2); 	% kg-m^2
lambda3 = (m/12)*(w^2 + L^2);  	% kg-m^2

J = diag([lambda1, lambda2, lambda3]);
%J = [2 -1 0; -1 4 0; 0 0 3];
orientation = SO3;
orientation.setToIdentity();
%w = 0.2*[1 2 2]';
w = [0.1, 20, 0.1]';
dt = 0.005;
T = orientation.m_data;
T(4,4) = 1.0;
%h = plottform(orientation.m_data);
view(135,30);
[h,pts] = scene.drawRobot3D(T,scale=scale);
scene.setBB(pts)
n = 1/dt;
for t = 0:dt:1
    wd = -inv(J)*(cross(w,J*w)); % angular acceleration by (3.14)
    w = w + wd * dt;
    orientation = orientation + w * dt; % update
    T = orientation.m_data';
    T(4,4) = 1.0;
    %plottform(T, handle=h);
    delete(h);
    h = scene.drawRobot3D(T,scale=scale);
    pause(3*dt);
end
