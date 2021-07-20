% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;
% given a plane
N1 = [1;1;1];
N1 = N1/norm(N1);
%given another plane
N2 = [0;0;1];
N2 = N2/norm(N2);
% given a Figure
O = [2,3,1; -2,3,1 ; -2,1,4; 2,1,4;];
O = O';
num = 30;
t = linspace(0,2*pi,num);
O = [cos(t);sin(t);zeros(1,num)];
O = rec3D(1,1,1);
O = [O;ones(1,size(O,2))];

% central orthogonal Transformation 
%T = (eye(3) - (1/dot(N1,N2))*N2*N1');
%Exemplo 2 
% trans = [1.5;1.5;1.5];
% T = translation(trans);
% O = T*O;
% T = translation(-trans);
% T = dilatation([1;1;-1])*T;
% TO = T*O;
% origin = [-0.5;-0.5;0.5];

%Exemplo 3
origin = [-0.5;-0.5;-0.5];
T = [0.877,0.479,0,0;-0.479,0.8772,0,0;0,0,1,0;0,0,0,1];
TO = T*O;





%% Template to make GIFS
filename = 'translationReflection.gif';
%Default configurations
fig = figure;
hold on

view(23,15);
%axis ([-1.5 2.5 -1.5 2.5 -1.5 2.5]);
axis equal
set(gca,'visible','off')
set(gcf,'color','white');
grid on
%camera
az = 23;
el = 15;
%% First frame
%drawPlan(N1,2.5);
%drawPlan(N2,4);
h = drawCube(O,'blue',0.05,0.3);
drawCoordinateSystem(origin  ,[1;0;0],[0;0;1]);
%plot3(O(1,:),O(2,:),O(3,:),'Color','blue');
frame = getframe(fig);
im{1} = frame2im(frame);
%% Intermediate Frames

im = Transformation(T,O,fig,im,az,el);

%% Last frame
%drawCoordinateSystem([0;0;0] - [1;1;0] ,[1;0;0],[0;0;1],2);
%plot3(TO(1,:),TO(2,:),TO(3,:),'Color','green');
drawCube(TO,'green');
frame = getframe(fig);
im{102} = frame2im(frame);

for idx = 1:102
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1.5);
    elseif (idx > 1) && (idx < 100)        
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.04);
    elseif idx == 102
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1.5);
    end

end



function im = Transformation (T,O,fig,im,az,el)
T3by3 = T(1:3,1:3);
T3by3(4,4) = 1;
Trans = inv(T3by3)*T;
nFrames = 101; % number of frames
% Método de Interpolação
[A,S,B] = svd(T(1:3,1:3));
U = B*S*B';
R = A*B';
U(4,4) = 1;
U = U*Trans;
eul = rotm2eul(R);
for i = 2:nFrames
    %R = eye(3) + (i/nFrames)*(T-eye(3));
    InStre = eye(4) + (i/nFrames)*(U-eye(4)); %Interpolate Stretch
    InAngle = rotm2eul(eye(3)) + (i/nFrames)*(eul-zeros(1,3));
    InRot = eul2rotm(InAngle);
    InRot(4,4) = 1;
    nv = InRot*InStre*O;
    %% Plot
    %h = fill3(nv(1,:),nv(2,:),nv(3,:),[0,1,0]);
    %h = plot3(nv(1,:),nv(2,:),nv(3,:),'Color','green');
    h = drawCube(nv,'green');
    %camera
    %view(az+i/3,el);
    drawnow
    frame = getframe(fig);
    im{i} = frame2im(frame);
    pause(0.01);
    delete(h);
    
    
end

end

function drawPlan(n,A,xo,color)
if nargin == 1
    A = 1;
    xo = [0;0;0];
    color = [1,0,0];
elseif nargin == 2
    xo = [0;0;0];
    color = [1,0,0];
end
[x,y] = findTriedro(n);
xp1 = xo + 2*A*x;
yp1 = xo + A*y;
xp2 = xo - 2*A*x;
yp2 = xo - A*y;
h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.2);
end
function [x,y] = findTriedro(z)
z = z/norm(z);
uTemp = [1;1;1];
uTemp = uTemp/norm(uTemp);
if ( 1 - dot(z,uTemp)  < 0.01)
    uTemp = [-uTemp(2);uTemp(1);uTemp(3)];
end
x = cross(uTemp,z);
y = cross(z,x);
end
function drawVector(x,color,p)
if nargin == 1
    p = [0;0;0];
    color = [0.2,0.4,0.8];
end
if nargin == 2
    p = [0;0;0];
end
quiver3(p(1),p(2),p(3),...
    x(1),x(2),x(3),'Color',color);
end
function O = rec3D(L,w,h)
unitCube = [-0.5 -0.5 -0.5; 0.5 -0.5 -0.5; -0.5 0.5 -0.5;0.5 0.5 -0.5; ...
    -0.5 -0.5 0.5;0.5 -0.5 0.5;0.5 0.5 0.5;-0.5 0.5 0.5];
a = -0.5;
b = 0.5;
n = 1;
x = a:(b-a)/n:b;
y = x;
z = x;
c = b*ones(1,n+1);
fac1 = [[-x;c;c] [-c;c;-z] [x;c;-c] ... 
        [c;c;z]];
fac2 = [[c;-y;c] [c;c;-z] [c;y;-c] ... 
        [c;-c;z]];
fac3 = [[-x;-c;c] [-c;-c;-z] [x;-c;-c] ... 
        [c;-c;z]];
fac4 = [[-c;-y;c] [-c;c;-z] [-c;y;-c] ... 
        [-c;-c;z]];
fac5 = [[-x;-c;c] [-c;-y;c] [x;c;c] ... 
        [c;y;c]];
fac6 = [[-x;-c;-c] [-c;-y;-c] [x;c;-c] ... 
        [c;y;-c]];
unitCube = [fac1 fac2 fac3 fac4 fac5 fac6];
I = eye(3);
T = I;
T(1:4:end) = [L w h];
O = T*unitCube;
end
function h = drawCube(O,color,faceAlpha,edgeAlpha)
if size(O,1) == 4
    O = O(1:3,:);
end
if nargin == 1
    color = [0 0 1];
    faceAlpha = 1;
    edgeAlpha = 1;
end
if nargin == 2
    faceAlpha = 1;
    edgeAlpha = 1;
end
[~,n]=size(O);
v = 1:(n/6); %6 faces
i = (n/6)*ones(1,(n/6));
fac = [v;v+i;v+2*i;v+3*i;v+4*i;v+5*i];
h = patch('Vertices',O','Faces',fac,...
    'FaceVertexCData',hsv(6),'FaceColor',color,'facealpha',faceAlpha,'edgealpha',edgeAlpha);
end
function M = translation(t)
M = eye(4);
M(1:3,4) = t;
end
function drawCoordinateSystem(p,x,z)
p = p - [0.1;0.1;0];
if nargin == 3
    index = [];
end
y = cross(z,x);
posTextX = p + x;
posTextY = p + y;
posTextZ = p + z;
text(posTextX(1),posTextX(2),posTextX(3),strcat('x',num2str(index)));
text(posTextY(1),posTextY(2),posTextY(3),strcat('y',num2str(index)));
text(posTextZ(1),posTextZ(2),posTextZ(3),strcat('z',num2str(index)));
quiver3(p(1),p(2),p(3),...
    x(1),x(2),x(3),'Color','red');
quiver3(p(1),p(2),p(3),...
    y(1),y(2),y(3),'Color','green');
quiver3(p(1),p(2),p(3),...
    z(1),z(2),z(3),'Color','blue');
text(p(1)-0.1*x(1),p(2),p(3),strcat('O',num2str(index)),'HorizontalAlignment','right');

end
function M = dilatation(s)
M = eye(4);
M(1,1) = s(1);
M(2,2) = s(2);
M(3,3) = s(3);
end
