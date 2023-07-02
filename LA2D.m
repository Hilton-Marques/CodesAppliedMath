clc;
clear all;
close all;


%Base do subespaço
b1 = [1;0];
b2 = [0;1];

O = [1.5, 0,-1.5,-1.5,0,1.5 ; ...
    sqrt(3)*[0.5,1,0.5, -0.5, -1, -0.5]];
O = Group(b1,b2);
c = min(O(2,:));

O = O + [0;-c];
O = box();
eixos = getBB(O);
figure
hold on
axis(eixos);
axis('off');
DrawObject(O);
%T = [2,0;0,2];
T = Translation([3,3])*Rotation(pi/1.5)*Shear(pi/4);
T = 1*quad2H(box()',trape()');
Transformation2(T,O,true)
function DrawObject(O,h)
h2 = fill(O(1,:), O(2,:), [0 0 1]);
set(h2, 'facealpha',0)
if (nargin > 1)
    h2 = fill(O(1,:), O(2,:), [0 0 1]);
    set(h2, 'facealpha',h)
end
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
end
function eixos =  getBB(O)
margin = 6.0;
x_max = max(O(1,:)) + margin;
x_min = min(O(1,:)) - margin;
y_max = max(O(2,:)) + margin;
y_min = min(O(2,:)) - margin;
eixos = [x_min,x_max,y_min,y_max];
end
function X = Transformation2(T,O,flag)
if nargin == 2
    flag = false;
end
T_old = T;
[~,c] = size(O);
nFrames = 100; % number of frames
% Método de Interpolação
X=O;
O(3,:) = 1.0;
X_end = T*O;
trans = [T(1,3),T(2,3)];
T(1,3) = 0.0;
T(2,3) = 0.0;
[A,S,B] = svd(T);
U = B*S*B';
R = A*B';
h2=fill(X(1,:), X(2,:), [0 0 1]);
%axis([-10 10 -10 10]);
set(h2, 'facealpha',0)
%pause(3);
for i = 1:nFrames
    fac = (i/nFrames);
    if flag
        Ti =  eye(3) + fac*(T_old-eye(3));
        nv = Ti*O;
    else        
        %Scale part
        Si = eye(3) + fac*(U-eye(3));
        %Rotation part
        eul = rotm2eul(T); %change rotation to euler angles
        InAngle = rotm2eul(eye(3)) + fac*(eul-zeros(1,3)); %interpo
        Ri = eul2rotm(InAngle); %bring back to matrix
        nv = Translation(fac*trans)*Ri*Si*O;
    end
    h2 = fill(nv(1,:)./nv(3,:), nv(2,:)./nv(3,:),[0 0 1]);
    set(h2, 'facealpha',0)
    pause(0.005);
    drawnow
    delete(h2);
end
h=fill(X_end(1,:)./X_end(3,:), X_end(2,:)./X_end(3,:), [0 0 1]);
set(h, 'facealpha',0);
end
function res = Rotation(teta)
res = eye(3);
res(1:2,1:2) = [[cos(teta),-sin(teta)];[sin(teta),cos(teta)]];
end
function res = Scale(fac)
res = eye(3);
res(1,1) = fac(1);
res(2,2) = fac(2);
end
function res = Translation(fac)
res = eye(3);
res(1,3) = fac(1);
res(2,3) = fac(2);
end
function M = quad2H(p_1,p_2)
n = [1,1,1];
z = n/norm(n);
y = cross(z,-[1,0,0]);
y = y/norm(y);
x = cross(y,z);
T = [x',y',n'];
p_1(:,3) = 1.0;
p_2(:,3) = 1.0;
H_1 = buildH(T*p_1');
H_2 = buildH(T*p_2');
M = inv(T)*H_2*inv(H_1)*T;
end
function H = buildH(p)
p0 = p(:,1);
p1 = p(:,2);
p2 = p(:,3);
p3 = p(:,4);
H = [p0,p1,p2];
lam = H\p3;
H = [lam(1)*p0,lam(2)*p1,lam(3)*p2];
end
function O = box()
O = ones(3,4);
O(1:2,1) = [-1,-1];
O(1:2,2) = [1,-1];
O(1:2,3) = [1,1];
O(1:2,4) = [-1,1];
end
function O = trape()
O = ones(3,4);
O(1:2,1) = [-4,-3];
O(1:2,2) = [4,-3];
O(1:2,3) = [1,2];
O(1:2,4) = [-2,2];
O = Translation([2.5,2.5])*Rotation(pi/3)*O;
end
function res = Shear(teta)
res = eye(3);
res(1:2,2) = [cos(teta),sin(teta)];
end
 