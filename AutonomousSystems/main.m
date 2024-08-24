%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Sunday, May 26th 2024, 3:03:25 pm
%Author: Hilton-Marques
%
%Description: This is a program to test the autonomous system provided in Farrell's
Aided Navigation book, example 3.16, pg. 90
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

clear all
clc;
close all;

addpath(genpath('../../my_libs'))

% Matrices
F = [[0.8187, 0.0000];[0.0906, 1.000]];
H = [0, 1];
G = [0.0906; 0.004];
L = [0.1832; 0.0187];

% Initial conditions
x0 = [10; 1.0];
e0 = x0;

% Create the autonomous system
sys = AutonomousSystem(F, H, L, G, x0, e0);
sys.Solver();
sys.exportFrame(filename="autonomous_linear");