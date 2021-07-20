clc;
clear all;
close all;


%Base do subespaço
b1 = [1;0];
b2 = [0;1];

%Subespaço
k = 10; % tamanho do subespaço
s = 100; % Número de linhas no grid
%grid (s,b1,b2,k);
%vector2DBasis(b1,[0 0 1])
%vector2DBasis(b2,[1 0 0])

% % Mapping non linear
% [D] = Object(); %Source points and Domain points
% TransformJ(D)
% F = [D(1,:).*cos(D(2,:));...
%     D(1,:).*sin(D(2,:))];
% % Input
% x =[2;3];
%Generate a orbit
T = Translation([3;0]);
S = Translation([1.5;(1.5*sqrt(3))]);
O = [1.5, 0,-1.5,-1.5,0,1.5 ; ...
    sqrt(3)*[0.5,1,0.5, -0.5, -1, -0.5]; ones(1,6)];
%Orbit(O,T,S,k);


O = Group(b1,b2); %Objeto a ser desenhado na configuração inicial

% Transformação
R = (sqrt(2)/2)*[1,-1,0;1,1,0;0,0,1];
O = R*O;
hold on
T = Translation([1;3]);
T = [1 0.5 0; 0.5 1 0;0 0 1];
M = R*T*R';
Transformation2(M,O);
keyboard;
O = T*O;
hold on
R = Reflection([],[1;-4]); 
%Q = Rotation(200*sqrt(2)*(2*pi)/360);
%T = Translation([1;2]);
Tf = Transformation2(R,O);

% for k = 1:20
%  Q = Rotation(k*100*sqrt(2)*(pi/180));
%  Transformation2(Q,O);
% end
% Xin = Transformation2(inv(T),X);


function space(x,b,k)
x = x/norm(x);
t = -k:k:k;
w = x(1)*t+b(1);
z = x(2)*t+b(2);
plot(w,z,'Color',[0.5 0.5 0.5],'Linewidth',1);
xlim([-k/2 k/2]);
ylim([-k/2 k/2]);
hold on
end

function gridBasis(l1,b1,b2,k)

for i = -l1:1:l1
    b1 = b1/norm(b1);
    b2 = b2/norm(b2);
    p_b1 = [-b1(2);b1(1)];
    p_b2 = [-b2(2);b2(1)];
    p_b1 = i*p_b1;
    p_b2 = i*p_b2;
    space(b1,p_b1,k);
    space(b2,p_b2,k);
end
end
function Thandle = triangle(xf,c,xi)
h = 0.25;
b = 0.25;
x = xf - xi;
p_x = [-x(2,1); x(1,1)];
p_x = p_x/norm(p_x);
p1 = xf + b/2*p_x;
p2 = xf - b/2*p_x;
p3 = xf + h*(x/norm(x));

X = [p1(1,1), p2(1,1) , p3(1,1)];
Y = [p1(2,1), p2(2,1) , p3(2,1)];
Thandle = fill(X, Y, c);
end

function vector2D(xf,c,xi)
if isempty(xf)
    line([0 xi(1,1)],[0 xi(2,1)], 'Color', [1 0 1],'Linewidth',2);
    triangle(xi,[1 0 1],[0;0]);
    return
end
if (nargin < 3)
    xi = zeros(2,1);
end
cx = xf(1,1)/norm(xf);
cy = xf(2,1)/norm(xf);
X = [xi(1,1), xf(1,1)];
Y = [xi(2,1), xf(2,1)];
line( X, Y, 'Color', c,'Linewidth',2);
line([0 xi(1,1)],[0 xi(2,1)], 'Color', [1 0 1],'Linewidth',2);
triangle(xf,c,xi);
triangle(xi,[1 0 1],[0;0]);
end

function vector2DBasis(x,c)
x = x/norm(x);
cx = x(1,1);
cy = x(2,1);

X = [0, x(1,1)];
Y = [0, x(2,1)];
line( X, Y, 'Color', c,'Linewidth',2);
triangle(x,c);
end

function M = Group(b1,b2)
c=0.2;
d = 1.5;
M = [ -1  (-1 + c)  (-1 +2*c)  (-1 +3*c)  (-1 +4*c)  (-1 +5*c)  (-1 +6*c)  (-1 +7*c) ...
    (-1 +8*c)  (-1 +9*c)  (-1 +10*c) (-1 +10*c) 1.5  1 1 -1 ; 2  (2-c) 2 (2-c) 2 (2-c) 2  (2-c)  2 (2-c)...
    2  (2 - d)  -1  -1 -2 -2];
M = [b1/norm(b1),b2/norm(b2)]*M; % Tem que se mostrar como seria nosso desenho na base canonica
%h=fill(M(1,:), M(2,:), [0 0 1]);
%set(h, 'facealpha',0)
[~,d] = size(M);
M = [M;ones(1,d)];
end

function X = Transformation2(T,O)
[~,c] = size(O);
nFrames = 100; % number of frames
% Método de Interpolação
%Render frame a frame
X=O;
h2=fill(X(1,:), X(2,:), [0 0 1]);
axis([-10 10 -10 10]);
set(h2, 'facealpha',0)
pause(3);
for i = 1:nFrames
    %Scale part
    Si = eye(3) + (i/nFrames)*(T-eye(3));
    %Rotation part
    %     eul = rotm2eul(T); %change rotation to euler angles
    %     InAngle = rotm2eul(eye(3)) + (i/nFrames)*(eul-zeros(1,3)); %interpo
    %     Ri = eul2rotm(InAngle); %bring back to matrix
    nv = Si*O;
    h2 = fill(nv(1,:), nv(2,:), [0 0 1]);
    set(h2, 'facealpha',0)
    pause(0.005);
    drawnow
    delete(h2);
end
X=T*O;
h=fill(X(1,:), X(2,:), [0 0 1]);
set(h, 'facealpha',0)
end

function DrawObject(O,h)

h2 = fill(O(1,:), O(2,:), [0 0 1]);
set(h2, 'facealpha',0)
if (nargin > 1)
    h2 = fill(O(1,:), O(2,:), [0 0 1]);
    set(h2, 'facealpha',h)
end


end

function [D] = Object()
o1 = 0.5;
o2 = 2.5;
r1 = 1;
r2 = 2;
n = 100;
y = o1:(o2-o1)/n:o2;
x = r1:(r2 -r1)/n:r2;
[X,Y] = meshgrid(x,y);
X = reshape(X',1,[]);
Y = reshape(Y',1,[]);
D = [X;Y]; % Domain points
%h2 = fill([a b],[a b], [0 0 1]);
%set(h2, 'facealpha',0)
end

%Transformation for Jacobian
function NGrid = TransformJ(D)
nFrames = 100; % number of frames
[~,n] = size(D);
D = [D(1,:);D(2,:);zeros(1,n)];
plot3(D(1,:),D(2,:),D(3,:),'o-');
%axis([(min(D(1,:))-1) (max(D(1,:))+1) (min(D(2,:))-1) (max(D(2,:))+1)]);
hold on
pause(1)
% Método de Interpolação
for i = 1:nFrames
    T = [(3*cos(D(1,:)) + cos(D(1,:)).*cos(D(2,:)));...
        (3*sin(D(1,:)) + sin(D(1,:)).*cos(D(2,:))); ...
        sin(D(2,:))];
    R = D + (i/nFrames)*(T-D);
    h2 = plot3(R(1,:),R(2,:),R(3,:));
    pause(0.01)
    delete(h2);
end
plot3(R(1,:),R(2,:),R(3,:));
end

function Q = Rotation(teta,center)
if nargin > 1
    T = Translation(center);
else
    center = [0;0];
end
Q = [cos(teta) -sin(teta) 0 ; sin(teta) cos(teta) 0; 0 0 1]; % Rotação
Q = T*Q;
end
function S = Scale(c)
S = [c 0 0;0 c 0; 0 0 1]; % Scale
end
function T = Translation(t)
T = [ [ 1 0 0]' [0 1 0]' [t;1]];
end
function R = Reflection(dir,normal)
%first argument is for line through the origin
%second argument is for lines whose normal is given
if nargin > 1
    dir = [-normal(2,1); normal(1,1)];
    dir = dir/vecnorm(dir);
    u = dir; 
    T = Translation(2*normal);
else
    u = dir;
    normal = 0*[-u(2,1); u(1,1)];
    T = eye(3);
end
%Plot line reflection
scaleFactor = 10;
xf = normal + scaleFactor*u;
xi = normal - scaleFactor*u;
grid on
line([xi(1,1) xf(1,1)],[xi(2,1) xf(2,1)],'Color','blue','LineStyle','--');
vector2D([],[0.5 0.5 0.5],normal);
hold on
I = eye(3);
I(3,3) = -1;
u = [u;0];
u_norm = vecnorm(u);
u = u/u_norm;
R = 2*(kron(u',u)) - I;
R = T*R;

end
function Orbit(O,T,S,k)
DrawObject(O);
for j = -k:k
    for i = -k:k
        newO = (T^i)*(S^j)*O;
        DrawObject(newO);
        if (i == 0) && (j == 0)
            DrawObject(newO,0.5);
        end
        hold on
        axis([-10 10  -10 10]);
    end
end
end
function Circle(r,center)
if nargin < 1
    center = [0;0];
end
t = 0:0.1:2*pi;
c = r*exp(complex(0,t));
c = c + complex(center); 
plot(c);
end
