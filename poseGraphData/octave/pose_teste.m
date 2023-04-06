clc;
clear all;
close all;

load intel-2d-posegraph.mat pg
disp(pg)
resErrorVec = edgeResidualErrors(pg);
plot(resErrorVec);