clear all
clc
close all

addpath("../my_libs/Manifolds/")

pose_a = SE2([0;0; pi/3]);

pose_b = SE2([4;4; 6*pi/3]);

Geodesic(pose_a, pose_b);
