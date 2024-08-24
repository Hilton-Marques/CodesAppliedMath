%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Tuesday, May 28th 2024, 7:21:10 pm
%Author: Hilton-Marques
%
%Description:A class to study naive and robust linear regression
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

close all
clear all
clc;

addpath(genpath("../../my_libs"))

%% Test Robust IEKF
obj = RobustEkf(method='l2', stdx=3.5, x_c=3.0, h=@(x) 0.27*(x - 3) + 0.27,H=0.27);

obj.ShowDistributionsPosterior(2.0, 1.0, false,false);

%% Test Iterated Weight
% Load the data
b = [1;2];

% Define the model 
H = [1; 1];

S = pinv(H);
x = S*b;
disp(x);

obj = IRWLS(H, b, "huber");
obj.Show(x);