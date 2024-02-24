clc;
clear;
close all;
rng('default');
addpath(genpath("../../my_libs/"));
%I am trying to reproduce this tweet:
% https://twitter.com/gabrielpeyre/status/1314068000912285696

% Export test images
obj = Mean();
obj.ShowCosts();
obj.exportFrame("filename","robust_distributions", "extension", "pdf");
% Gaussian sums
obj.GaussianSum();
obj.exportFrame("filename","gaussian_sum", "extension", "pdf");

% Generate data
n = 100;
s = 1;
x = s*rand(1, n) - s/2;
y = s*rand(1, n) - s/2;

pts = [x; 0*y];
m_prev = mean(pts,2);
var_prev = var(pts(1,:));
pts_close = pts;
%obj.ShowDataGaussian(pts_close, obj.m_blue);

n_outliers = 10;
t = 5*[1;1];
s = 10;
x = s*rand(1, n_outliers) - s/2 + t(1);
y = s*rand(1, n_outliers) - s/2 + t(2);
pts_outliers = [x; 0*y];
%obj.ShowDataGaussian(pts_outliers, obj.m_red);
pts = [pts, pts_outliers];

%% Empirical determination of std increasing
m = mean(pts,2);
delta = m - m_prev; %increasing in mean
n_tot = (n_outliers + n);
eps = n_outliers/n_tot;
std = (delta(1)/eps);
var_new = (n/n_tot)*var_prev + (n_outliers/n_tot)*std^2;
var_ref = var(pts(1,:));
%% 

%pts = [14, 4; 9, 5; 11, 7]';
obj = Mean(pts);
obj.ShowCosts();
obj.ShowPts(pts_close, obj.m_red);
obj.ShowPts(pts_outliers, obj.m_blue);


%obj.exportFrame(filename='robust_inference', extension='pdf');
%obj.Solver();
%x = obj.m_x;
%Visual data
x = [0.0732; -0.0225];
x = [0.1160; 0];
plot(m(1), m(2), 'o', 'MarkerFaceColor', obj.m_green, 'markersize', 6, 'MarkerEdgeColor','black');
axis("tight");
%obj.exportFrame(filename='l2_norm', extension='pdf');
plot(x(1), x(2), 'o', 'MarkerFaceColor', 'cyan', 'MarkerEdgeColor','black', 'markersize', 6);
obj.exportFrame(filename='robust_inference', extension='pdf');


