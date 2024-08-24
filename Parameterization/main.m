%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Wednesday, April 17th 2024, 2:23:12 pm
%Author: Hilton-Marques
%
%Description: Parameterizations of S2
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}


clc;
clear all;
close all;
addpath(genpath("../../my_libs/"));

obj = GeodesicCircle();
P = [1,0,0];
P = P/norm(P);
Q = [0,0,1];
obj.Solver(P,Q);