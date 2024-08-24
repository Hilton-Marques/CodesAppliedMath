%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Wednesday, March 13th 2024, 6:30:55 pm
%Author: Hilton-Marques
%
%Description: Here we test the Bode plot of the transfer function in Maybeck (1979), page 302.
This transfer function arises from a simple error-state kalman filter problem.
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

%Clear configuration
clear all;
close all;
clc;

%Transfer function
s = tf('s');
wn = 2;
G = -sqrt(2) * wn * (s + wn/sqrt(2)) / (s*s + sqrt(2)* wn * s + wn*wn);
%bode(G);

plotoptions = bodeoptions;
plotoptions.Grid = 'on';
plotoptions.FreqScale = 'linear';
plotoptions.Title.String = 'Bode Plot of Transfer Function';
bodeplot(G,plotoptions,{0, 5 * wn});
exportgraphics(gcf,'bode_plot_ekf.jpeg');
