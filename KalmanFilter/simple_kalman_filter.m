%{
%Copyright (c) 2024 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Thursday, May 9th 2024, 10:30:07 am
%Author: Hilton-Marques
%
%Description:
%HISTORY: Simple example of 1d kalman filter
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}
close all;
clc;
clear all;

figure;
hold on;
%define two real symbolic variables
syms x y

%Prediction
%x_c = 20; %barfoot h3
%std_x = 3; %barfoot h3
x_c = 3.0; %h2
std_x = 3.5; %h2
p_x = (1/(std_x * sqrt(2*pi))) * exp(((x - x_c).^2 ./ (std_x^2)) * -0.5);
e1 = double(p_mean(x, p_x));

%Measurement
%std_y = 0.3; % barfoot h3
%yk = 0.938461538461539; %barfoot h3

y_c = 2.0; 
std_y = 1;
p_yx = (1/(std_y * sqrt(2*pi))) * exp(((y - h2(x))^2 ./ (std_y ^ 2)) * -0.5);
p = p_x * p_yx;

%Exact Posterior
p_post = subs(p, y, y_c);
p_norm = int(p_post, x, -inf,inf);
p_post = p_post / p_norm;

e2 = double(p_mean(x, p_post));

% Calculate the mode
% Derive the first derivative with respect to y
dp_dx = diff(p_post, x);

% Solve for zeros of the derivative to find critical points
critical_points = solve(dp_dx == 0, x);

% Convert critical points to numeric if they are not already
numeric_critical_points = double(critical_points);

%Measurement Linearized
p_yx_lin = (1/(std_y * sqrt(2*pi))) * exp(((y - h2_lin(x,x_c))^2 ./ (std_y ^ 2)) * -0.5);
p_lin = p_x * p_yx_lin;

% Linearized Posterior
p_post_lin = subs(p_lin, y, y_c);
p_norm_lin = int(p_post_lin, x, -inf,inf);
p_post_lin = p_post_lin / p_norm_lin;

e3 = double(p_mean(x, p_post_lin));

%Apply IEKF
[v,ps,p_joint] = iekf(e3, y, std_y, p_x, y_c, x);

e4 = v(end);
% Find the closest critical point to e4
[a,b] = min(abs(numeric_critical_points - e4));
% The MAP estimation is 
map = numeric_critical_points(b);
error = abs(map - e4);
fprintf('Compare MAP with IEKF %d', error)

if (true)
	x_inv = linspace(map - 8, map + 4);
	fac = 10;
	hh = [];
	for i = 1:height(v)
		color = rand(3,1);
		xline(v(i), '--','LineWidth',1)
		hh = [plot(x_inv, fac*subs(ps(i), x, x_inv) + y_c, 'color', color, 'linewidth',1.5), hh];
		show_contour(map, y_c, 8, p_joint(i), 0.02,color);
		%hh = [plot(x_inv, h2_lin(x_inv, v(i)),'-','color',color,'linewidth',1.5),hh];
		%hh = [plot(v(i), h2(v(i)), 'o', 'MarkerFaceColor','black','markersize',8)];
	end
end

% Plot the contour plot
red = [0.9176 0.2627 0.2078];
blue = [0.2588 0.5216 0.9569];
green = "#699C52";


%title('Contour plot of p(x, y)');
%xlabel('x');
%ylabel('y');
x_inv = linspace(map - 8, map + 4);
p1 = plot(x_inv, h2(x_inv),'-','color','#CF5044','linewidth',1.5);
%p2 = plot(x_inv, h2_lin(x_inv, x_c),'-','color','#1C758A','linewidth',1.5);
fac = 10;
%p3 = plot(x_inv, fac*subs(p_post, x, x_inv) + y_c,'color',red,'linewidth',1.5);
p4 = plot(x_inv, fac*subs(p_post_lin, x, x_inv) + y_c, 'color', blue, 'linewidth',1.5);
p5 = plot(x_inv, fac*subs(p_x, x, x_inv) - 2*y_c,'color',green,'LineWidth',1.5);

%plot contours
%show_contour(map, y_c, 8, p, 0.02,red);
show_contour(map, y_c, 8, p_lin, 0.02, blue);
yline(y_c, '--','LineWidth',1)
xline(map,'--','LineWidth',1)
xline(e3,'--','LineWidth',1) 
%xline(e2,'--','LineWidth',1) %mean
xline(x_c,'--','LineWidth',1) %prior

%legend([p1, p2, p5,  p3, p4], '$y=h(x)$','$y=\bar{h}(x)$','$p_0(x|y)$','$p(x|y)$','$\bar{p}(x|y)$','interpreter','latex','location','northwest');
legend([p1, p5, p4, hh(1), hh(2), hh(3)], '$y=h(x)$', '$p_0(x|y)$', '$\bar{p}_0(x|y)$','$\bar{p}_1(x|y)$','$\bar{p}_2(x|y)$','$\bar{p}_3(x|y)$','interpreter','latex','location','northwest');

xlabel('$x$ (estado)','interpreter','latex');
ylabel('$y$ (medicao)','interpreter','latex');
exportgraphics(gcf,'iterative_iekf.pdf','ContentType','vector')
exportgraphics(gcf,'iterative_iekf.png','Resolution','300')


keyboard

function y = h(x)
	y = x.^2/20;
end

function y = h2(x)
	y = 0.01 * x.^ 3; 
end

function e = p_mean(x, p)
e = int(x * p, x, -inf, inf);
end

function y = derh2(x)
	y = 0.01 * 3 * x .^ 2;
end

function y = h2_lin(x, xo)
	y = h2(xo) + derh2(xo) * (x - xo);
end

%From barfoot, pg. 97
function y = h3(x)
f = 400;
b = 0.1;
y = f*b/x;
end

function y = derh3(x)
f = 400;
b = 0.1;
y = f*b * -1 * x.^-2;
end

function y = h3_lin(x, xo)
	y = h3(xo) + derh3(xo) * (x - xo);
end

function [v,ps, p_joint] = iekf(xo, y, std_y, p_x, yk, x, flag_plot)
max_iter = 3;
v = zeros(max_iter,1);
ps = {};
p_joint = [];

for i = 1:max_iter
	v(i) = xo;
	p_yx_lin = (1/(std_y * sqrt(2*pi))) * exp(((y - h2_lin(x,xo))^2 ./ (std_y ^ 2)) * -0.5);
	p_lin = p_x * p_yx_lin;
	% Linearized Posterior
	p_post_lin = subs(p_lin, y, yk);
	p_norm_lin = int(p_post_lin, x, -inf,inf);
	p_post_lin = p_post_lin / p_norm_lin;
	ps = [ps, p_post_lin];
	p_joint = [p_joint, p_lin];
	xo = double(p_mean(x, p_post_lin));
end

end

function show_contour(x_c, y_c, margin, p, z_value,color)
n = 300;
x_inv = linspace(x_c - margin, x_c + 4, n);
y_inv = linspace(y_c - margin, x_c + 4, n);
[x_grid, y_grid] = meshgrid(x_inv,y_inv);
p_func = matlabFunction(p);
z = p_func(x_grid, y_grid);
contourf(x_grid, y_grid, z, [z_value, z_value],'FaceAlpha',0.5,'FaceColor',color); % 50 contour levels
end