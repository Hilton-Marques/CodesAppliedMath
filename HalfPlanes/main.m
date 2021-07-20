clc;
clear;
close all;


p0 = [0.0,0.0,0.0] - [0.5,0.5,0.5] ;
p1 = [1.0,0.0,0.0] - [0.5,0.5,0.5] ;
p2 = [1.0,1.0,0.0] - [0.5,0.5,0.5] ;
p3 = [0.0,1.0,0.0] - [0.5,0.5,0.5] ;

p4 = [0.0,0.0,1.0] - [0.5,0.5,0.5] ;
p5 = [0.0,1.0,1.0] - [0.5,0.5,0.5];
p6 = [1.0,1.0,1.0] - [0.5,0.5,0.5];
p7 = [1.0,0.0,1.0] - [0.5,0.5,0.5];

v0 = Vertex(p0',1);
v1 = Vertex(p1',2);
v2 = Vertex(p2',3);
v3 = Vertex(p3',4);
v4 = Vertex(p4',5);
v5 = Vertex(p5',6);
v6 = Vertex(p6',7);
v7 = Vertex(p7',8);
hf1 = HalfPlane([v0;v1;v2;v3]);
hf2 = HalfPlane([v4;v5;v6;v7]);
hf3 = HalfPlane([v2;v1;v7;v6]);
hf4 = HalfPlane([v3;v5;v4;v0]);
hf5 = HalfPlane([v4;v7;v1;v0]);
hf6 = HalfPlane([v2;v6;v5;v3]);
hfs1 = [hf1,hf2,hf3,hf4,hf5,hf6];
coords1 = [p0;p1;p2;p3;p4;p5;p6;p7];
coords2 = zeros(8,3);
pts1 = [v0,v1,v2,v3,v4,v5,v6,v7];
pts2(length(pts1)) = Vertex();
hfs2(length(hfs1)) = HalfPlane();
teta = pi/3;
trans = 1*[0.5,0.5,0.5];
m = [cos(teta),sin(teta),0,trans(1);-sin(teta),cos(teta),0,trans(2);0,0,1,trans(3);...
    0,0,0,1];
for i = 1:length(pts1)
    coord1 = pts1(i).coord;
    coord1(4,1) = 1.0;
    newCoord = m*coord1;
    coords2(i,:) = newCoord(1:3,1)';
    pts2(i) = Vertex(newCoord(1:3,1),i);
end
for i = 1:length(hfs1)
    hfi = hfs1(i);
    hfs2(i) = HalfPlane([pts2(hfi.pts(1).id);pts2(hfi.pts(2).id);...
        pts2(hfi.pts(3).id);pts2(hfi.pts(4).id)]);
end
solid1 = Solid(hfs1,pts1);
solid2 = Solid(hfs2,pts2);
%start gif
filename = 'interseção entre cubos.gif';
%Default configurations
fig = figure;
hold on
angleInit = 30;
view(angleInit,30);
axis equal;
set(gca,'visible','off');
set(gcf,'color','white');

solid1.show('b');
solid2.show('r');

pause(0.01);
count = 1;
frame = getframe(fig);
im{count} = frame2im(frame);
%[fRot,count,im] = rot(angleInit,count,im,fig);
%find chebshive point
% 
% A = [0.0124932555  , 0.0401622131  , 0.999115050  ;...
% -0.0124932565 , -0.0401622131 , -0.999115050 ;...
% -0.882352948  , 0.470588237 , -0.00000000    ;...
% 0.857492924  , -0.514495790 , 0.00000000      ;...
% 0.498283893  , .867013931 , -0.00000000      ;...
% -0.417733222  , -0.908569753 , 0.00000000    ;...		
% -0.0545788631 , 0.00602463400 , 0.998491287 ;...
% 0.0545788631  , -0.00602463400 , -0.998491287;...
% -0.482935488  , 0.875655949 , -0.00000000    ;...
% 0.529282928  , -0.848445475 , 0.00000000      ;...
% 0.813164353  , .582034051 , -0.00000000      ;...
% -0.863017797  , -0.505173624 , 0.00000000    ];
% 
% A(:,4) = 1.0;
% 
% b = -[2976.4672851562500;...
% -2983.4611816406250;...
% -735.44116210937500;...
% 692.42553710937500;...
% -1130.9051513671875;...
% 1094.3043212890625;...
% 3009.2534179687500;...
% -3016.2429199218750;...
% -1387.2362060546875;...
% 1262.4964599609375;...
% -650.15728759765625;...
% 433.45388793945312];

[A,b] = getAb(solid1,solid2);
%A = A(1:6,:);
%b = b(1:6,:);
n = size(A,2);
m = size(A,1);
c = zeros(n,1);
c(end) = -1;
tic
center = linprog(c,A,b);
r = center(4,1);
center = center(1:3);
count = 0;
figure (2)
hold on
view(30,30);

plot3(center(1),center(2),center(3),'o','MarkerSize',5,'MarkerFaceColor','r');

for i = 1:12
    n = A(i,1:3);
    d = b(i);
    
    pos = -n*d;
    u = center' - pos;
    plot3(pos(1),pos(2),pos(3),'o');
    v = n;
    %drawPlan(n',d,'r');
    %quiver3(pos(1),pos(2),pos(3),1000*n(1),1000*n(2),1000*n(3));
    x = u / norm(u);
   % quiver3(pos(1),pos(2),pos(3),u(1),u(2),u(3));
    value = dot(x,n);
    if (value > 0)
        count = count + 1;
    end
end
% for i = 6:12
%     n = A(i,1:3);
%     d = b(i);
%     v = n;
%     drawPlan(n',d,'b');
%     u = center' + d*n ;
%     value = dot(u,n);
%     if (value < 0)
%         count = count + 1;
%     end
% end
% [X,Y,Z] = sphere;
% X = r*X + center(1);
% Y = r*Y + center(2);
% Z = r*Z + center(3);
% plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','r');
% [fRot,count,im] = rot(fRot,count,im,fig);
% surf(X,Y,Z)
% [fRot,count,im] = rot(fRot,count,im,fig);
% createGif(im,filename);
% keyboard
newA = A;
ids = find (b < 0);
negIds = find(b>=0);
nIds = length(ids);
idsM = zeros(size(A,1),nIds);
for i = 1:nIds
    idsM(ids(i),i) = -1.0;
end
newA(ids,:) = -newA(ids,:);
cAux = zeros(m,1);
I = eye(m,m);
newff = [0*c;0*-c;zeros(nIds,1);cAux(negIds);ones(nIds,1)];
newAA = [newA,-newA,idsM,I(1:m,negIds),-idsM];
newA = [newA,-newA,idsM,I];
cAux(ids) = 1.0;
newf = [0*c;0*-c;zeros(nIds,1);cAux];
b(ids) = - b(ids);
n = size(newf,1);
basis = zeros(1,m);
count = 0;
count2 = 1;
for i = 1 : m 
    flag = false;
    for j =  1 : length(ids)
        if i == ids(j)
            basis(i) = 20 + j+count;
            ids(j) = [];
            flag = true;
            count = count + 1;
            break;
        end
    end
    if flag
        continue;
    end
    basis(i) = count2 + 8 + nIds;
    count2 = count2 + 1;
end
tic
%[x,points,cost,basis] = solveLP(newA,b,newf,[n-12+1:n],b);
[x,points,cost,basis] = solveLP(newAA,b,newff,basis,b);
[xOther,fval] = linprog(newff,-eye(n,n),zeros(n,1),newAA,b);
toc
%basis = sort(basis);
init = x(basis);
newf = zeros(20,1);
n = size(newf,1);
newAA(:,end-nIds+1:end) = [];
newf(4,1) = -1.0; 
newf(8,1) = 1.0; 
tic
[x,points,cost,basis] = solveLP(newAA,b,newf,basis,init);
toc

%plot3(x(1,1),x(2,1),x(3,1),'o','MarkerSize',10,'MarkerFaceColor','r');
% tic
% [xOther,fval] = linprog(newf,-eye(n,n),zeros(n,1),newA,b([negIds;ids]));
% toc
% newA = [A,-A,eye(m,m)];
% newf = [c;-c;ones(m,1)];
% newN = size(newA,2);
% tic
% x = linprog(newf,-eye(newN,newN),zeros(newN,1),newA,b);
% toc
% center = x(1:3,1);
% 
% plot3(x(1),x(2),x(3),'o','MarkerSize',10,'MarkerFaceColor','r');

pause(0.01);
%[fRot,count,im] = rot(fRot,count,im,fig);

radius = x(end);
solid1.translateHfs(center);
solid2.translateHfs(center);
pointsForFecho = getPointsForFecho(solid1,solid2);
h2 = showPoints(pointsForFecho);

pause(0.01);
%[fRot,count,im] = rot(fRot,count,im,fig);

[k1,av1] = convhull(pointsForFecho);
ids = unique(k1(:));
h = trisurf(k1,pointsForFecho(:,1),pointsForFecho(:,2),pointsForFecho(:,3),'FaceColor','cyan');
pause(0.01);
%[fRot,count,im] = rot(fRot,count,im,fig);

dualPoints = zeros(size(k1,1),3);
triAreas = [];
dds = [];
totalArea = 0;
for i = 1:size(k1,1)
    inc = k1(i,:);
    p1 = pointsForFecho(inc(1),:);
    p2 = pointsForFecho(inc(2),:);
    p3 = pointsForFecho(inc(3),:);
    normal = cross(p2 - p1, p3 - p1);
    normN  = norm(normal);
    normal = normal / norm(normal);
    d = dot(normal,p1);
    dualPoints(i,:) = (1/d)*normal + center';
    dds(end+1) = d;
    totalArea = totalArea + normN/2;
end
ptsFecho = pointsForFecho(ids,:);
normals = [];
hfsSolidIntersection(size(ptsFecho,1)) = HalfPlane();

for i = 1:size(ptsFecho,1)
    p = ptsFecho(i,:);
    triAreas(end+1) = 1/norm(p);
    normals(end+1,1:3) = p/norm(p);
    hfsSolidIntersection(i) = HalfPlane([],normals(i),triAreas(i));
end
triAreas = sort(triAreas);
%showPoints(dualPoints);
newAreas = [];
[k2,av1] = convhull(dualPoints);
dds2 = [];
normalsTri = [];
for i = 1:size(k2,1)
    inc = k2(i,:);
    p1 = dualPoints(inc(1),:);
    p2 = dualPoints(inc(2),:);
    p3 = dualPoints(inc(3),:);
    normal = cross(p2 - p1, p3 - p1);
    normN = norm(normal);
    normal = normal / normN;
    d = dot(normal,p1);
    dds2(end+1) = 1/d;
    normalsTri(end + 1, 1:3) = normal;
    %newAreas(end+1) = normN/2;
    newAreas(end+1) = d;
end
areasPlan = zeros(6,1);
areasPlan(1) = newAreas(1) + newAreas(3);
areasPlan(2) = newAreas(4) + newAreas(5);
areasPlan(3) = newAreas(6) + newAreas(9);
areasPlan(4) = newAreas(11) + newAreas(12);
areasPlan(5) = newAreas(7) + newAreas(10);
areasPlan(6) = newAreas(2) + newAreas(8);
newAreas = sort(newAreas);
delete(h2);
delete(h);
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','cyan');
pause(0.01);
[fRot,count,im] = rot(fRot,count,im,fig);

function [A,b] = getAb(solid1,solid2)
[A1,b1] = solid1.getA();
[A2,b2] = solid2.getA();
A = [A1;A2];
b = [b1;b2];
end
function pts = getPointsForFecho(solid1,solid2)
pts1 = solid1.getPointsForFecho();
pts2 = solid2.getPointsForFecho();
pts = [pts1;pts2];
end
function h2 = showPoints(points)
h2 = plot3(points(1,1),points(1,2),points(1,3),'o','color','black');
for i = 1:size(points,1)
    p = points(i,:);
    h2(end+1) = plot3(p(1),p(2),p(3),'o','color','black');
    
end
end
function [fRot,count,im] = rot(initangle,count,im,fig)
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
for i = 1:2:180
    fRot = initangle + i;
    view(fRot,30);
    pause(0.05);
    count = count + 1;
    frame = getframe(fig);
    im{count} = frame2im(frame);
end
end

function [x,points,cost,basis] = solveLP(A,b,c,basis,initPoint)
points = [];
At = A';
x = zeros(size(c,1),1);
d = x;
dAux = x;
B = A(:,basis);
Binv = inv(B);
x(basis) = initPoint; % xpos at starting corner
cost = c(basis)'*x(basis);     % cost at starting corner
n = size(A,2);

points(end+1,1:3) = x(1:3,1) - x(4:6,1);
for iter = 1:100
    d = zeros(size(c,1),1);
    %dAux = zeros(size(c,1),1);
    y = Binv' * c(basis);           % this y may not be feasible
    nonBasics = setdiff(1:n,basis);
    for in = nonBasics
        rmin = c(in) - At(in,:)*y;
        if rmin < -.00000001
            break
        end
    end
    %[rmin,in] = min(c - A'*y); % minimum r and its index in
    if rmin >= -.00000001      % optimality is reached, r>=0
        check = A*x - b;
        break;                 % current x and y are optimal
    end
    Aj = A(:,in);
    d(basis) = Binv * Aj;  % decrease in x from 1 unit of xin
    dAux(basis) =  -d(basis);
    dAux(in) = 1;
    A*dAux
    xb = x(basis);
    db = d(basis);
    ids = find(db > 0 );
    tetas = xb(ids)./db(ids);
    [teta,~] = min(tetas);
    %Blands rule
    if teta == 0
        outs = find(teta==xb);
    else
        outs = ids(find(teta==tetas));
    end
    out = outs(1);
%     for i = 2:size(outs)
%         outTemp = outs(i);
%         if basis(outTemp) < basis(out)
%             out = outTemp;
%         end
%     end
    cost = cost + teta*rmin;  % lower cost at end of step
    x(basis) = x(basis) - teta*d(basis);   % update old x
    x(in) = teta;      % find new positive component of x
    check = A*x - b
    basis(out) = in;      % replace old index by new in basis
    B = A(:,basis);
    BinvL = (1/db(out))*Binv(out,:);
    Q = zeros(12,12);
    BinvBef = Binv;
    for i = 1:12
        if ( i == out)
            Binv(i,:) = (1/db(i))*Binv(i,:);
            %Q(i,i) =  (1/db(i));
        else
            Binv(i,:) = Binv(i,:) + (-db(i))*BinvL;
            Q(i,i) =  1.0;
            Q(i,out) =  ( -db(i)/db(out) );
        end
    end
    err = norm(Binv - inv(B));
    points(end+1,1:3) = x(1:3,1) - x(4:6,1);
    Binv = inv(B);
end

end
function createGif(im,filename)
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
    elseif (idx > 1) && (idx < length(im))
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
    end
end
end

function drawPlan(n,d,color,A,xo)
if nargin == 1
    A = 2000;
    xo = [0;0;0];
    color = [1,0,0];
    d = 0;
elseif nargin == 2
    xo = [0;0;0];
    color = [1,0,0];
    A = 2000;
elseif nargin == 3
    xo = [0;0;0];
    A = 2000;
end
[x,y] = findTriedro(n);
xp1 = xo + A*x - n*d;
yp1 = xo + A*y - n*d;
xp2 = xo - A*x - n*d;
yp2 = xo - A*y - n*d;
h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.5);
end
function [x,y] = findTriedro(z)
z = z/norm(z);
uTemp = [z(3);z(1);-z(2)];
uTemp = uTemp/norm(uTemp);
if (dot(z,uTemp) == 1)
    uTemp = uTemp([2,1,3]);
end
x = cross(uTemp,z);
y = cross(z,x);
end