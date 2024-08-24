%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Sunday, April 28th 2024, 10:59:59 am
%Author: Hilton-Marques
%
%Description:
%HISTORY: Testing the adjoint matrix of the SE2 group
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}
clear all;
clc;
close all;

addpath(genpath("../../my_libs/"))

t = [1;2;3];
v = [4;5;6];
x = SO3();
K = x + v;
-inv(K.m_data)*SO3.hat(t)
-SO3.hat(inv(K.m_data)*t)*inv(K.m_data)

x = SO3();
v = [-0.0002644277,0.0002056362,-0.0003052109];
a = x + v;
a.m_data
Y = SE2([3,3,pi/2]);
u = 0.5*[1;1;0.5];

adj = Adjoint();

%adj.ConjugateMap(u, Y);
%adj.LieBracketApprox(u, Y);
A = SE2([4,-0.8,0]);
B = SE2([3,2, pi/9]);
adj.FrameTransformation(A,B);
%adj.adjoint_rep(A,B);

