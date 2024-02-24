%{
%Copyright (c) 2023 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Sunday, November 12th 2023, 11:48:43 pm
%Author: Hilton-Marques
%
%Description: A educational program to understand continuous transformations
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}
clc;
clear all;
close all;
addpath(genpath("../../my_libs/"));
Matrices();
%LinearTransformations()

function RotationGroup()
T = @(X,Y) {-Y, X};
obj = ContinuousTransformations("rotation_group.gif");
obj.IntegrateInfinitesimalTransformation(T);
obj.save(repeat=true);
end

function Matrices()
M = 1/2*[[5,1];[1,5]];
L = logm(M);
T = @(X,Y) {L(1,1) * X + L(1,2) * Y, L(2,1) * X + L(2,2) * Y};
obj = ContinuousTransformations("symmetric_matrices_polar.gif");
obj.IntegrateInfinitesimalTransformation(T);
obj.save(repeat=true);
end

function LinearTransformations()
M = 1/2*[[5,1];[1,5]];
L = (M - 3*eye(2));
T = @(X,Y) {L(1,1) * X + L(1,2) * Y, L(2,1) * X + L(2,2) * Y};
obj = ContinuousTransformations("linear_matrices.gif");
obj.LinearTransformation(T);
obj.save(repeat=true);
end

function R = logmapp(A, n)
M = A - eye(2);
R = zeros(2);
for i = 1:n
    R = R + (-1)^(i+1) * (M)^i / i;
end
end