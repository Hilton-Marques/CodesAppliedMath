clear all
close all
clc;
rng('default');

addpath(genpath("../my_libs"));

d1 = [1;0;0];
d2 = [1;0;0];
d3 = [1;0;0];

obj = Solver(d1,d2,d3,[-1;1;3]);
obj.setTarget([1;-1;3]);
obj.Compute();
obj.setTarget([1;-1;-3]);
obj.Compute();
obj.save(true);