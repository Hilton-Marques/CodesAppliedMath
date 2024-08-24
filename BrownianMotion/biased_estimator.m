%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Saturday, February 24th 2024, 11:05:14 pm
%Author: Hilton-Marques
%
%Description: Biased estimator using a random walk as bias
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

clc;
clear all;
close all;


%Parameters
n = 500; %Number of steps
x = 0; %Initial position
z_std = 1; %Standard deviation of the noise
z_mean = 0; %Mean of the noise
m = 1; %realizations
y = zeros(m, n); %sample path of the process or realization
x_std = 1.5;


%Generate the sample path
y(:,1) = x;
t = 1:n;
for i = 1:m
    for j = 1:n-1
        y(i, j+1) = normrnd(0, x_std) +  (y(i, j) + normrnd(z_mean,z_std));
    end
end

%Plot the sample path
figure
hold on
plot(repmat(t,m,1)', y', '-o','MarkerFaceColor',"auto",'LineWidth',1.5,'markersize',3,'MarkerEdgeColor','none');
h= plot(repmat(t,m,1), y, 'o','MarkerFaceColor',"auto","markersize",3, 'MarkerEdgeColor', 'black');
%set(h, {'MarkerFaceColor'}, get(h,'Color'));
%set(h, {'MarkerEdgeColor'}, 'black');

means = mean(y')
mean(y)
figure
plot(t, mean(y), 'color', 'black', 'linewidth',2)
exportgraphics(gcf,'GaussianWalk.png','Resolution',300);