%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Sunday, March 10th 2024, 11:29:54 am
%Author: Hilton-Marques
%
%Description: A program to test signal that are preserved by LTI systems
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

% Clear the workspace
clear all;
close all;
clc;

% Add the path to Scene
addpath("../../../Projetos/my_libs/");

y_1 = Signal(r=1.0,omega=1.0,phase=pi/3);
y_2 = Signal(r=1.01, omega=-1.0, amplitude=1.5);
%y_3 = Signal(0.8,-2,1.5);

signals = Signals();
signals = signals + y_1;
%signals = signals + y_1 + y_2;
%signals.Show();
signals.ShowDynamically();


