clc;
clear;
close all;

addpath(genpath("../my_libs/"))
obj = Pendulum();

%propose some function for theta
theta_0 = pi/3;
theta_fcn = @(t) theta_0 - 2 * theta_0*t;
%theta_fcn = @(t) theta_0*cos(t * pi);
obj.Trajectory(theta_fcn);
obj.save();