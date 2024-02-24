clear all
clc
close all

addpath("../../my_libs/Manifolds/")

pose_a = SE3([0;0;15;1;1;1]);

pose_b = SE3([0;0;1;0;0;0]);

rot_a = SO3(eye(3));
rot_b = SO3(pose_a.m_data(1:3, 1:3));

Geodesic(rot_a, rot_b);
