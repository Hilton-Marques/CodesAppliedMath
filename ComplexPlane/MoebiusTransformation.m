clc;
clear all;
close all;

[X,Y,Z] = Sphere();
SphereObj = [X Y Z];
x = -1:0.2:1;
[X,Y]=meshgrid(x);
Of = complex(X,Y);
angle = [0,0,-pi/3];
R = eul2rotm(angle); R(4,4) = 1;
R(1:3,4) = [0;0;2];
%R = [1 0 0 0;0 1 0 0;0 0 1 3; 0 0 0 1];
Obj1 = {SphereObj Of};
Obj2 = InvProject(Obj1);
Obj3 = TransformSphere(R,Obj2);
Obj4 = Project(Obj3);
Obj1 = {SphereObj [X Y 0*X]};
Objs = {Obj1 Obj2 Obj3 Obj4};
set(gca,'visible','off')
h = DrawObject(Obj1);
delete(h(1));
for i = 1:3
  Render(Objs{i},Objs{i+1});
end
DrawObject(Obj4);

function FObj = InvProject(O)
%The Inverse Projection is done on the initial sphere
SphereObj = O{1};
[n,~] = size(SphereObj);
northPole = [SphereObj(1,1),SphereObj(1,(n+1)),SphereObj(1,(2*n+1))];
%Objects
s = length(O);
FObj{1} = SphereObj;
IObj{1} = SphereObj;
for i = 2:s
P = O{i};
[n,~] = size(P);
X = real(P);
Y = imag(P);
Z = zeros(n,n);
P = [X Y Z];
vecPoints = [X(:) Y(:) Z(:)];
xo = vecPoints;
xf = repmat(northPole,n^2,1);
u = xf - xo;
m = vecnorm(u')';
%x= 2*sin(teta)
s = 2*(2./m);
factor = m - s;
spherePoints = xo + factor.*(u./m);
FObj{i} = [reshape(spherePoints(:,1),n,n)...
    reshape(spherePoints(:,2),n,n)...
    reshape(spherePoints(:,3),n,n)];
IObj{i} = P;
end
%Render(IObj,FObj);
end
function Sterio(plan)
%Initialize sphere
teta = (0:pi/12:2*pi)';
n = length(teta);
fi = (0:pi/(n-1):pi)';
X = sin(fi)*(cos(teta))';
Y = sin(fi)*(sin(teta))';
Z = ((cos(fi))) + 1;
Z = repmat(Z,1,n);
surf(X,Y,Z)
%Find points that verify plan equation X - Y = 0
p = plan(X,Y,Z);
c1 = find(abs(p) < 0.01);
nc = size(c1,1);
Xlevel = X(c1);
Ylevel = Y(c1);
Zlevel = Z(c1);
points = [Xlevel Ylevel Zlevel];
northPole = [0 0 2];
%Build lines from North Pole to points
scaleFactor = 3;
linex = [((-0.5*(points(:,1) - northPole(1)) ) + northPole(1)) ...
    ((scaleFactor*(points(:,1) - northPole(1)) ) + northPole(1))];
liney = [((-0.5*(points(:,2) - northPole(2)) ) + northPole(2)) ...
    ((scaleFactor*(points(:,2) - northPole(2)) ) + northPole(2))];
linez = [((-0.5*(points(:,3) - northPole(3)) ) + northPole(3)) ...
    ((scaleFactor*(points(:,3) - northPole(3)) ) + northPole(3))];
%Find points from line that go to complex plane
%director vector
xo = repmat(northPole,nc,1);
v = points - xo;
%proportion rule
scale = 2 * 1 ./ (2*ones(nc,1) - points(:,3));
xf = xo + (scale .* v);
%avoid infinits or zeros
c1 = find(xf(:,3) == 0);
c2 = find((points(:,3) ~= 0 & points(:,3) ~= 2));
plot3(xf(c1,1), xf(c1,2), xf(c1,3),'Color','red','Marker','o');
plot3(points(c2,1),points(c2,2), points(c2,3),'Color','yellow','Marker','o');
%plot lines from north pole to points
for i = 1:nc
    line(linex(i,:),liney(i,:),linez(i,:));
end
axis([-3 3 -3 3 0 2.5]);
drawnow
end
function Render(I,F)
%clf;
s = length(I);
nFrames = 20;
hold on
for i = 1:s
    P = I{i};
    [n,~] = size(P);
    h(i) = surf(P(:,1:n),P(:,(n+1):(2*n)),P(:,2*n+1:end));
    if (i == 1)
        set(h,'FaceAlpha',0.5,'LineStyle','none');
    end
    axis([-8 3 -8 3 0 5]);
    view(19,28);
end
delete(h);
for i = 1:nFrames
    for j = 1:s
        T = F{j};
        P = I{j};
        [n,~] = size(T);
        Ti = P + (i/nFrames)*(T - P);
        h(j) = surf(Ti(:,1:n),Ti(:,(n+1):(2*n)),Ti(:,2*n+1:end));
        if (j == 1)
            set(h(j),'FaceAlpha',0.5,'LineStyle','none');
        end
    end
    drawnow;
    pause(0.1);
    delete(h);
end
for i = 1:s
    T = F{i};
    [n,~] = size(T);
    h(i) = surf(T(:,1:n),T(:,(n+1):(2*n)),T(:,2*n+1:end));
    if (i == 1)
        set(h(i),'FaceAlpha',0.5,'LineStyle','none');
    end
end
delete(h);
hold off
end
function [X,Y,Z] = Sphere()
teta = (0:pi/12:2*pi)';
n = length(teta);
fi = (0:pi/(n-1):pi)';
X = sin(fi)*(cos(teta))';
Y = sin(fi)*(sin(teta))';
Z = ((cos(fi))) + 1;
Z = repmat(Z,1,n);
end
function FObj = TransformSphere(R,O)
s = length(O);
for i = 1:s
    P = O{i};
    [n,~] = size(P);
    %Take the coordinate points
    X = P(:,1:n);
    Y = P(:,(n+1):(2*n));
    Z = P(:,2*n+1:end);
    Obj = [X(:)';Y(:)';Z(:)';ones(1,n^2)];
    Obj = R*Obj;
    FObj{i} = [reshape(Obj(1,:),n,n) ...
              reshape(Obj(2,:),n,n) ...
              reshape(Obj(3,:),n,n)];
end
%Render(O,FObj);
end
function FObj = Project(O)
s = length(O);
SphereObj = O{1};
FObj{1} = SphereObj;
[n,~] = size(SphereObj);
%The northPole is always the first element of the matrix
northPole = [SphereObj(1,1),SphereObj(1,(n+1)),SphereObj(1,(2*n+1))];
for i = 2:s
    points = O{i};
    [n,~] = size(points);
    X = points(:,1:n);
    Y = points(:,(n+1):(2*n));
    Z = points(:,2*n+1:end);
    vecPoints = [X(:) Y(:) Z(:)];
    %director vector
    xo = repmat(northPole,n^2,1);
    v = vecPoints - xo;
    %proportion rule
    scale = northPole(3) ./ (northPole(3) - vecPoints(:,3));
    xf = xo + (scale .* v);
    FObj{i} = [reshape(xf(:,1),n,n) reshape(xf(:,2),n,n) reshape(xf(:,3),n,n)];
    
end
%Render(O,FObj);
end
function h = DrawObject(O)
%clf;
s = length(O);
hold on
for i = 1:s
    P = O{i};
    [n,~] = size(P);
    h(i) = surf(P(:,1:n),P(:,(n+1):(2*n)),P(:,2*n+1:end),'FaceColor','interp');
    if (i == 1)
        set(h(i),'FaceAlpha',0.5,'LineStyle','none');
    end
    axis([-8 3 -8 3 0 5]);
    view(19,28);
end
pause(1)
end
