% Clear workspace
clear
close(findall(0,'Type','figure'));
clc
% Parameters

n  = 3;    
ids = [1,2,5];
m  = 1000*rand(1,n); 
m = [500,550,600,650,600];
m_new = m(ids);
x  = (rand(3,n)-2*rand(3,n));     
x = [[-1,1,-1,1];[1,-1,-1,1];[1,1,-1,-1]];
x(:,end+1) =  [3;3;3];
x_mean = mean(x(:,1:end-1),2);
x_new = [x(:,1),x_mean, x(:,5)];
v  = zeros(3,n); 
dt = .0008; 
colors = [[0.7963 0.4416 0.4462];...
    [0.4657 0.2790 0.6754];...
    [0.9037 0.9085 0.7472];...
    [0.2605 0.6896 0.1318]];
colors(end+1,:) = [0.2588,0.5216,0.9569];
colors_new = colors(ids,:);
x = x_new;
m = m_new;
colors = colors_new;
% Graphics
% 
% fig = figure( ...
%     'color', [.5, .5, .7], ...
%     'menubar', 'none', ...
%     'numbertitle', 'off', ...
%     'name', 'n-body');                             
% 
% ax = axes( ...
%     'clipping','off', 'box', 'on',      'ticklength',[0,0], ...
%     'xticklabel',[],  'yticklabel',[],  'zticklabel',[], ...
%     'xlim', [-4,4],   'ylim', [-4,4],   'zlim', [-4,4], ...
%     'xgrid', 'on',    'ygrid', 'on',    'zgrid', 'on', ...
%     'color',[0,0,0],  'gridcolor',[1,1,0]); 
gif_obj = Gif(strcat('Nbody','.gif'));
gif_obj.setAngle(1);
hold on
for e = 1:n                                            
    planets(e) = line( ...
         'marker', 'o', ...
         'markersize', m(e)/50, ...
         'markerfacecolor',colors(e,:));
end
expansion(planets,x,gif_obj);
scene(planets,x,gif_obj);

% Loop

tic;

while toc < 100
    %view(10*toc,30);   
    %gif_obj.update()
    planets = plotUpdate(planets,x);       
    gif_obj.update()
    [x,v]   = dataUpdate(x,v,m,dt);                    
    %pause(0.1); 
end
gif_obj.save();
% Plot Update Function

function planets = plotUpdate(planets,x)       
    hold on
    for e = 1:length(x)
        set(planets(e), ...
            'XData',x(1,e), ...
            'YData',x(2,e), ...
            'ZData',x(3,e));
        plot3(x(1,e),x(2,e),x(3,e),'o','markerfacecolor','blue','markersize',0.5);
    end
end

% Data Update Function 

function [x,v] = dataUpdate(x,v,m,dt)                   
    for i = 1:length(x)                                   
        a(:,i) = zeros(3,1);
        for j = 1:length(m)                                 
            if i ~= j
                a(:,i) = a(:,i) + ...
                    m(j).*(x(:,j)-x(:,i))./(norm(x(:,j)-x(:,i))^3+.5);
            end
            v(:,i) = v(:,i) + a(:,i)*dt;
            x(:,i) = x(:,i) + v(:,i)*dt;
        end
    end
end

function h = vec(p1,p2,color,fac)
L = norm(p2 - p1);
axis = p2 - p1;
p_trans = p1 + fac*axis;
h = line([p1(1),p_trans(1)],[p1(2),p_trans(2)],...
                [p1(3),p_trans(3)],'Color',color,'Linewidth',1);

axis = axis / norm(axis);
axis_rot = cross(axis,[0,0,1]);
angle =  acos(dot(axis,[0,0,1]));
rot_matrix = axang2rotm([axis_rot,-angle]);


hold on
[X,Y,Z]=cylinder([0 (1)],20 );
[m,n] = size(X);
p = -0.08*[X(:)';Y(:)';Z(:)'];
p = rot_matrix * p;
X = reshape(p(1,:),m,n) + p_trans(1);
Y = reshape(p(2,:),m,n) + p_trans(2);
Z = reshape(p(3,:),m,n) + p_trans(3);
h(end+1) = surf(X,Y,Z,'FaceColor','black','EdgeColor','blue');
end
function scene(planets,x,gifobj)
out = bb(x);
axis(out);
plotUpdate(planets,x)  
gifobj.exportFrame('inital_pose');
n = 4;
xs = x(:,end);
h = [];
for i = 1:4
    xf = x(:,i);
    h_i = vec(xs,xf,'cyan',0.9);
    h(:,end+length(h_i)) = h_i;
end
gifobj.exportFrame('interation');
delete(h)
c = mean(x(:,1:4),2);
hold on
plot3(c(1),c(2),c(3),'o','markersize',18,'markerfacecolor','magenta');
h = [];
for i = 1:4
    xf = x(:,i);
    h_i = vec(xf,c,[0.9176,0.2627,0.2078],0.7);
    h(:,end+length(h_i)) = h_i;
end
gifobj.exportFrame('pole');
delete(h);
vec(xs,c,'cyan',0.9);
gifobj.exportFrame('final_1');
plot3(3.5,3.5,-1.5,'o','markersize',10,'markerfacecolor','blue');
vec([3.5;3.5;-1.5],c,'cyan',0.9);
gifobj.exportFrame('final_2');

end
function expansion(planets,x,gifobj)
hold on
out = bb(x,4.5);
axis(out);
plotUpdate(planets,x)  
vec(x(:,3),x(:,2),'cyan',0.9);
vec(x(:,1),x(:,2),'cyan',0.9);
text(x(1,1),x(2,1),x(3,1),'x_f','color','white','VerticalAlignment','top');
text(x(1,2),x(2,2),x(3,2),'x_c','color','white','VerticalAlignment','top');
text(x(1,3),x(2,3),x(3,3),'x_s','color','white','VerticalAlignment','top');
[X,Y,Z] = sphere(20);
r = norm(x(:,3) - x(:,2));
X2 = X * r ;
Y2 = Y * r;
Z2 = Z * r;
surf(X2,Y2,Z2,'EdgeAlpha',0.1,'FaceColor','white','FaceAlpha',0.2,'EdgeColor','white');
R = Rot3D(pi,[0,0,1]');
vec(x(:,2),R*x(:,3),'red',0.9);
u = 0.5*(R*x(:,3) - x(:,2));
text(u(1),u(2),u(3),'{\it r}','fontsize',15,'color','white')
gifobj.exportFrame('expansion_r');
end
function out = bb(x,margin)
xmin = min(x(1,1:end)) - margin;
xmax = max(x(1,1:end)) + margin;
ymin = min(x(2,1:end)) - margin;
ymax = max(x(2,1:end)) + margin;
zmin = min(x(3,1:end)) - margin;
zmax = max(x(3,1:end)) + margin;
out = [xmin,xmax,ymin,ymax,zmin,zmax];
end
function Rot3D = Rot3D(angle, axis)
n = axis / vecnorm(axis); 
chute = [n(3,1);-n(1,1);n(2,1)];
plan1 = (eye(3) - n*n')*chute;
plan1 = plan1/ vecnorm(plan1); 
plan2 = cross(plan1,n);
R = [cos(angle) -sin(angle) 0;...
    sin(angle) cos(angle) 0; ...
    0 0 1];
M = [plan1 plan2 n];
Rot3D = M*R*M';
end