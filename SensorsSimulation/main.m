clear all;
close all;
clc;

addpath(genpath("../../../Projetos/RVC3-MATLAB/toolbox"));
addpath(genpath("../../../Projetos/my_libs"));
rng('default')
[TrueMotion,imu] = imudata();
l = 0.301;    % body length
w = 0.088;    % body width
b = 0.05;     % body height
dt = imu.dt;
wt = imu.gyro;
obj = IMU(l,w,b,wt,dt,imu.accel, imu.magno);
%Show bias
obj.ShowBias();


%Covariance propgation
rots = obj.Solver(true, TrueMotion.orientation);
%obj.CompareWithExact(rots,TrueMotion.orientation,"error_raw");
%rots2 = obj.SolverCMF(true, TrueMotion);
%obj.CompareWithExact(rots2, TrueMotion.orientation,"error_cmp");
keyboard;