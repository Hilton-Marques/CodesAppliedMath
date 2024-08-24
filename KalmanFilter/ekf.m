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
std_x = 1; %h2
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

%Measurement Linearized
p_yx_lin = (1/(std_y * sqrt(2*pi))) * exp(((y - h2_lin(x,x_c))^2 ./ (std_y ^ 2)) * -0.5);
p_lin = p_x * p_yx_lin;

% Linearized Posterior
p_post_lin = subs(p_lin, y, y_c);
p_norm_lin = int(p_post_lin, x, -inf,inf);
p_post_lin = p_post_lin / p_norm_lin;

e3 = double(p_mean(x, p_post_lin));

% Plot the contour plot
red = [0.9176 0.2627 0.2078];
blue = [0.2588 0.5216 0.9569];
green = "#699C52";

%title('Contour plot of p(x, y)');
%xlabel('x');
%ylabel('y');
x_inv = linspace(e3 - 8, e3 + 4);
%p1 = plot(x_inv, h2(x_inv),'-','color','#CF5044','linewidth',1.5);
p2 = plot(x_inv, h2_lin(x_inv, x_c),'-','color','#1C758A','linewidth',1.5);
fac = 10;
%p3 = plot(x_inv, fac*subs(p_post, x, x_inv) + y_c,'color',red,'linewidth',1.5);
p4 = plot(x_inv, fac*subs(p_post_lin, x, x_inv) + y_c, 'color', blue, 'linewidth',1.5);
p5 = plot(x_inv, fac*subs(p_x, x, x_inv) - 2*y_c,'color',green,'LineWidth',1.5); %prior

%plot contours
%show_contour(map, y_c, 8, p, 0.02,red);
%show_contour(e3, y_c, 8, p_lin, 0.02, blue);
yline(y_c, '--','LineWidth',1)
xline(e3,'--','LineWidth',1)
yline(h2(x_c),'--','LineWidth',1)
xline(x_c,'--','LineWidth',1)
%legend([p1, p2, p5,  p3, p4], '$y=h(x)$','$y=\bar{h}(x)$','$p_0(x|y)$','$p(x|y)$','$\bar{p}(x|y)$','interpreter','latex','location','northwest');
%legend([p1, p5, p4, hh(1), hh(2), hh(3)], '$y=h(x)$', '$p_0(x|y)$', '$\bar{p}_0(x|y)$','$\bar{p}_1(x|y)$','$\bar{p}_2(x|y)$','$\bar{p}_3(x|y)$','interpreter','latex','location','northwest');

% xlabel('$x$ (estado)','interpreter','latex');
% ylabel('$y$ (medicao)','interpreter','latex');
% exportgraphics(gcf,'iterative_iekf.pdf','ContentType','vector')
% exportgraphics(gcf,'iterative_iekf.png','Resolution','300')

% Chech mahalanobis distance
c = [x_c; h2(x_c)];
xhat = [e3; y_c];
v = xhat - c;
z_value = eval_dist(p_lin, xhat);
%show_contour(e3, y_c, 8, p_lin, z_value, blue);

%tangent point
%plot(xhat(1), xhat(2), 'o','markersize',6,'MarkerFaceColor','black');
hh = derh2(x_c);
hh = [1;hh];
angle = atan2(hh(2), hh(1));
p0 = c + -10*hh ;
p1 = c + 10*hh;

%Compare with mahalanobis distance of 1
%show_contour(e3, y_c, 8, p_lin, 0.028, blue);
S = [[std_x^2, derh2(x_c)*std_x^2];[derh2(x_c)*std_x^2, std_y^2 + derh2(x_c)^2 * std_x^2]];
%ShowEllipse(inv(S),b,2.1*angle);

%Normal equation
S = [[std_x^2, 0];[0, std_y^2]];
%S = eye(2);
Cinv = inv(S);
L  = chol(S);
A_star = inv(L')*hh;
check = pinv(A_star);
H = hh;
b = [0;y_c - h2(x_c)];
P_i = H'*Cinv*H;
new_meas = P_i^-1*H'*Cinv*b;
new_Y = H*P_i^-1*H'*Cinv*b;
r = sqrt((new_Y - b)'*Cinv*(new_Y - b));
% Conditional Expectation?

ShowEllipse(inv(S),c + b, r);
plot(xhat(1), h2(x_c) + new_Y(2), 'o','markersize',6,'MarkerFaceColor','black');
plot(c(1), c(2), 'o','markersize',6,'MarkerFaceColor','black');
plot(x_c, -4, 'o','markersize',6,'MarkerFaceColor','black');
plot(c(1)+new_meas, y_c, 'o','markersize',6,'MarkerFaceColor','black');
arrow([c;0.4],[c + b;0.4],'color', '#699C52','tipWidth', 0.060,'stemWidth', 0.035);
arrow([c+b;0],[c + (new_Y);0],'color', '#644172','tipWidth', 0.060,'stemWidth', 0.035);
xlim([-2,10])
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
camlight
exportgraphics(gcf,"only_map.png","Resolution",300)
axis equal
%line([p0(1),p1(1)],[p0(2), p1(2)]);
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
z = contourf(x_grid, y_grid, z, [z_value, z_value],'FaceAlpha',0.5,'FaceColor',color); % 50 contour levels
end

function z = eval_dist(p, x)
p_func = matlabFunction(p);
z = p_func(x(1), x(2));
end

function ShowEllipse(A, t, r)
L = chol(A)';
T = inv(L);
n = 100;
theta = linspace(0,2*pi,n);
x = r*[cos(theta); sin(theta)];
%rot = [[cos(angle), -sin(angle)];[sin(angle), cos(angle)]];
y = T * x + t;
fill(y(1,:), y(2,:), 'red','FaceAlpha',0.5);
end

function h = arrow(p1,p2,varargin)
%mArrow3 - plot a 3D arrow as patch object (cylinder+cone)
%
% syntax:   h = mArrow3(p1,p2)
%           h = mArrow3(p1,p2,'propertyName',propertyValue,...)
%
% with:     p1:         starting point
%           p2:         end point
%           properties: 'color':      color according to MATLAB specification
%                                     (see MATLAB help item 'ColorSpec')
%                       'stemWidth':  width of the line
%                       'tipWidth':   width of the cone
%
%           Additionally, you can specify any patch object properties. (For
%           example, you can make the arrow semitransparent by using
%           'facealpha'.)
%
% example1: h = mArrow3([0 0 0],[1 1 1])
%           (Draws an arrow from [0 0 0] to [1 1 1] with default properties.)
%
% example2: h = mArrow3([0 0 0],[1 1 1],'color','red','stemWidth',0.02,'facealpha',0.5)
%           (Draws a red semitransparent arrow with a stem width of 0.02 units.)
%
% hint:     use light to achieve 3D impression
%
propertyNames = {'edgeColor'};
propertyValues = {'none'};
%% evaluate property specifications
for argno = 1:2:nargin-2
	switch varargin{argno}
		case 'color'
			propertyNames = {propertyNames{:},'facecolor'};
			propertyValues = {propertyValues{:},varargin{argno+1}};
		case 'stemWidth'
			if isreal(varargin{argno+1})
				stemWidth = varargin{argno+1};
			else
				warning('mArrow3:stemWidth','stemWidth must be a real number');
			end
		case 'tipWidth'
			if isreal(varargin{argno+1})
				tipWidth = varargin{argno+1};
			else
				warning('mArrow3:tipWidth','tipWidth must be a real number');
			end
		otherwise
			propertyNames = {propertyNames{:},varargin{argno}};
			propertyValues = {propertyValues{:},varargin{argno+1}};
	end
end
%% default parameters
if ~exist('stemWidth','var')
	ax = axis;
	if numel(ax)==4
		stemWidth = norm(ax([2 4])-ax([1 3]))/300;
	elseif numel(ax)==6
		stemWidth = norm(ax([2 4 6])-ax([1 3 5]))/300;
	end
end
if ~exist('tipWidth','var')
	tipWidth = 3*stemWidth;
end
tipAngle = 22.5/180*pi;
tipLength = tipWidth/tan(tipAngle/2);
ppsc = 50;  % (points per small circle)
ppbc = 250; % (points per big circle)
%% ensure column vectors
p1 = p1(:);
p2 = p2(:);
%% basic lengths and vectors
x = (p2-p1)/norm(p2-p1); % (unit vector in arrow direction)
y = cross(x,[0;0;1]);    % (y and z are unit vectors orthogonal to arrow)
if norm(y)<0.1
	y = cross(x,[0;1;0]);
end
y = y/norm(y);
z = cross(x,y);
z = z/norm(z);
%% basic angles
theta = 0:2*pi/ppsc:2*pi; % (list of angles from 0 to 2*pi for small circle)
sintheta = sin(theta);
costheta = cos(theta);
upsilon = 0:2*pi/ppbc:2*pi; % (list of angles from 0 to 2*pi for big circle)
sinupsilon = sin(upsilon);
cosupsilon = cos(upsilon);
%% initialize face matrix
f = NaN([ppsc+ppbc+2 ppbc+1]);
%% normal arrow
if norm(p2-p1)>tipLength
	% vertices of the first stem circle
	for idx = 1:ppsc+1
		v(idx,:) = p1 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
	end
	% vertices of the second stem circle
	p3 = p2-tipLength*x;
	for idx = 1:ppsc+1
		v(ppsc+1+idx,:) = p3 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
	end
	% vertices of the tip circle
	for idx = 1:ppbc+1
		v(2*ppsc+2+idx,:) = p3 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
	end
	% vertex of the tiptip
	v(2*ppsc+ppbc+4,:) = p2;
	% face of the stem circle
	f(1,1:ppsc+1) = 1:ppsc+1;
	% faces of the stem cylinder
	for idx = 1:ppsc
		f(1+idx,1:4) = [idx idx+1 ppsc+1+idx+1 ppsc+1+idx];
	end
	% face of the tip circle
	f(ppsc+2,:) = 2*ppsc+3:(2*ppsc+3)+ppbc;
	% faces of the tip cone
	for idx = 1:ppbc
		f(ppsc+2+idx,1:3) = [2*ppsc+2+idx 2*ppsc+2+idx+1 2*ppsc+ppbc+4];
	end
	%% only cone v
else
	tipWidth = 2*sin(tipAngle/2)*norm(p2-p1);
	% vertices of the tip circle
	for idx = 1:ppbc+1
		v(idx,:) = p1 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
	end
	% vertex of the tiptip
	v(ppbc+2,:) = p2;
	% face of the tip circle
	f(1,:) = 1:ppbc+1;
	% faces of the tip cone
	for idx = 1:ppbc
		f(1+idx,1:3) = [idx idx+1 ppbc+2];
	end
end
%% draw
fv.faces = f;
fv.vertices = v;
h = patch(fv);
for propno = 1:numel(propertyNames)
	try
		set(h,propertyNames{propno},propertyValues{propno});
	catch
		disp(lasterr)
	end
end
end
