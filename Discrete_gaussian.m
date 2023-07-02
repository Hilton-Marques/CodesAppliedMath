clc;
clear;
close all;
sigma = 1.25
x = linspace(-5,5,9);

g = exp((-x.^2)/(2*sigma^2));
g = g/sum(g)

a = fspecial('gaussian',[1,9],1.0)