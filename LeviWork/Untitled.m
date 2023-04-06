clc;
clear;
close all;

x = linspace(0,1,100);
y = f(x);



function res = f(x)
res = log(x ./ (1 - x));
end