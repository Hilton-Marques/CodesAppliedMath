%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Sunday, February 25th 2024, 9:58:25 pm
%Author: Hilton-Marques
%
%Description: We describe a brownian motion by sampling 
% its derivative by a N(0, ti - ti-1) distribution.
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

clc;
clear all;
close all;

dt = 0.01;
t = 1.0;
n = t/dt;
m = 10; % realizations
z_mean = 0;
z_std = sqrt(dt);
x = 0; %Initial position
y = zeros(m, n);

% Main loop
y(:,1) = x;
for i = 1:m
  for j = 2:n
    y(i, j) = y(i, j-1) + normrnd(z_mean, z_std);
  end
end


%Plot the sample path
t = linspace(0, t, n);

figure
hold on
plot(repmat(t,m,1)', y', '-o','MarkerFaceColor',"auto",'LineWidth',1.5,'markersize',3,'MarkerEdgeColor','none');
h= plot(repmat(t,m,1), y, 'o','MarkerFaceColor',"auto","markersize",3, 'MarkerEdgeColor', 'black');
set(h, {'MarkerFaceColor'}, get(h,'Color'));
%set(h, {'MarkerEdgeColor'}, 'black');
exportgraphics(gcf,'brownian_motion.png','Resolution',300);
figure
plot(t, mean(y), 'color', 'black', 'linewidth',2)
exportgraphics(gcf,'brownian_mean.png','Resolution',300);