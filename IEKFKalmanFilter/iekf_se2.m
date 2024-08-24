%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Monday, March 4th 2024, 6:23:28 pm
%Author: Hilton-Marques
%
%Description: Here we formulate the covariance propagating algorithm in the SE2 group from 
% https://github.com/UMich-CURLY-teaching/UMich-ROB-530-public/tree/main/code-examples/MATLAB/matrix_groups
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

clear all;
close all;
clc;

% Add the path to SE2 library
addpath(genpath("../../../Projetos/my_libs"));

% generate a path
n = 10;
x_init = SO3();
v = x_init.getRandomTangent();
v = 0.1 * ones(1,3);

%Initialize robots to be simulated
n_simu = 200;
x_simu(n_simu) = SO3();

% Init as elements of so2
x(n) = SO3();
for i = 1:n-1
    x(i+1) = x_init + v;
    x_init = x(i);
end

%Init control inputs (or measruments)
u = zeros(3, n);
for i = 1:n-1
  u(:,i) = x(i+1) - x(i);
end

%Noise parameters
%Q = diag([0.03^2, 0.03^2, 0.1^2]);
Q = 0.0001 * eye(3);
% Cholesky factor of covariance for sampling
L = chol(Q, 'lower');


green = [0.2980 .6 0];
figure
hold on
plot(0, 0,'.','Color', [green, .25], 'markersize', 14);
view(30,30);
%Define process noise
for i = 1:n-1
    x_simu = propagation(x_simu, u(:,i), L);
    p = zeros(2, n_simu);
    %p = zeros(3, n_simu);
    for j = 1:n_simu
        t = x_simu(j).m_data(1:2,3);
        t = rotm2axang(x_simu(j).m_data);
        p(1,j) = t(1);
        p(2,j) = t(2);
        p(3,j) = t(3);
    end
    % show particles
    %plot(p(1,:), p(2,:), 'o', 'Color',green);
    t = [i;0,;0];
    p = p + t;
    plot3(p(1, :),p(2, :),p(3, :), 'o', 'color', green);
end


function x_simu = propagation(x_simu, u, std)
  for i = 1:size(x_simu,2)
      x_simu(i) = ((x_simu(i) + u) + x_simu.getRandomLieAlgebra(std));
  end
end

