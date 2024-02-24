clear all
clc
close all

addpath("../../my_libs/Manifolds/")

pose_a = SE3([0;0;15;1;1;1]);

pose_b = SE3([0;0;1;0;0;0]);

Geodesic(pose_a, pose_b);
