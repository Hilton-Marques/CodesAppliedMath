clc;
clear;
close all;

%% Halfplanes
conec = [0, 3, 2; 0, 1, 3; 
         4, 1, 0; 4, 5, 1; 
         5, 3, 1; 5, 7, 3;...
         4, 2, 6; 4, 0, 2;
         7, 4, 6; 7, 5, 4;
         6, 3, 7; 6, 2, 3];
conec_fitting = [0,1,3,2;...
          4,6,7,5;...
          5,7,3,1;...
          4,0,2,6;...
          4,5,1,0;...
          2,3,7,6];
      
case0 = [0, 0, 10; 10, 0, 10; 0, 10, 10; 10, 10, 10;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
case1 = [0, 0, 10;0, 0, 10;0, 0, 10;0, 0, 10;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
case2 = [0, 0, 10;10, 0, 10;0, 10, 10;10, 10, 10;0, 0, 0;10, 0, 10;0, 10, 0;10, 10, 0];
case3 = [0, 0, 10;10, 0, 10;0, 10, 10;10, 10, 10;0, 0, 0;10, 0, 10;0, 10, 10;10, 10, 0];
case4 = [0, 0, 0;10, 0, 10;0, 10, 10;10, 10, 0;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
case5 = [0, 0, 0;0, 0, 0;0, 10, 10;10, 10, 0;0, 0, 0;0, 0, 0;0, 10, 0;10, 10, 0];




case_bruno =  [1276.25464 -2885.52686 -3043.70361;...
     1376.55835 -2885.95142 -3046.62451;...
     1279.13794 -2785.50562 -3049.44556;...
     1379.58325 -2786.01099 -3058.30469;...
     1276.25464 -2885.52686 -3043.70361;...
     1376.55835 -2885.95142 -3046.62451;...
     1279.13794 -2785.50562 -3049.44556;...
     1379.58325 -2786.01099 -3058.30469];
 pts = {case1,case2,case3,case4,case5};
 real_volumes = [];
 for i = 1:5     
     figure
     hold on
     view(-172,33);
     volume = calculateVolume(pts{i},conec+1,i);
     
     real_volumes(end + 1) = volume;
     name = strcat('caseDegeneratedVolume',num2str(i));
     title(name);
     exportgraphics(gca,strcat(name,'.png'),'Resolution',1000);
     hold off
 end
 for i = 1:5
    out = getVectorizedPoints(conec_fitting,pts{i});
 end


A = [-0.00000000 , -0.00000000 , 1.00000000, -10.0000000;...
0.00000000  , 0.00000000  , -1.00000000 , 0;...
0.707106829 , 0.00000000  , 0.707106829, -7.07106829;...
-1.00000000 , 0.00000000  , 0.00000000 ,  0;...
0.00000000  , -1.00000000 , 0.00000000 , 0;...
0.00000000  , 0.707106829 , 0.707106829, -7.07106829];

B = [0.00000000  , 0.00000000   , 1.00000000   , -10;...
0.367001027 , -0.367001027 , -0.854763448 , 2.13690853;...
1.00000000  , 0.00000000   , 0.00000000   , -10.0000000;...
-1.00000000 , 0.00000000   , 0.00000000   , 0;...
0.00000000  , -1.00000000  , 0.00000000   , 0;...
0.00000000  , 1.00000000   , 0.00000000   , -10.0000000];

C = [0.00000000  , 0.00000000  , 1.00000000 , -10.0000000;...
0.00000000  , 0.00000000  , -1.00000000, 5.00000000;...
1.00000000  , 0.00000000  , 0.00000000 , -10.0000000;...
-1.00000000 , 0.00000000  , 0.00000000 , 0;...
0.00000000  , -1.00000000 , 0.00000000 , 0;...
0.00000000  , 1.00000000  , 0.00000000 , -10.0000000];

D = [0.00000000  0.00000000  , 1.00000000  , -5.0;...
0.00000000  0.00000000  , -1.00000000 , 0;...
1.00000000  0.00000000  , 0.00000000  , -10;...
-1.00000000 0.00000000  , 0.00000000  , 0;...
0.00000000  -1.00000000 , 0.00000000  , 0;...
0.00000000  1.00000000  , 0.00000000  , -10];

E = [0.577350259 , -0.577350259 , 0.577350259, 0;...
0.00000000  , 0.00000000   , -1.00000000, 0;...
-1.00000000  , -0.00000000  , -0.00000000, 0;...
-1.00000000  , -0.00000000  , -0.00000000, 0;...
-0.00000000 , -1.00000000  , -0.00000000, 0;...
0.00000000  , 1.00000000   , 0.00000000 , -10];

normals = {A,B,C,D,E};

volumes = [];
for i = 1:5
    if i == 5
        a = 1;
    end
    pt_i = pts{i};

    figure
    hold on
    view(-172,33);    
    %axis ([-1,11,-1,11,-1,11]);
    plotSolid(pt_i,conec + ones(12,3),'cyan','blue');
    [k2,dualPoints,volume] = findInterior(normals{i},[12,0,0]);
    volumes(end+1) = volume;
    %plotPlane(normals{i});
    plot3(dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'o','MarkerFaceColor','yellow');
    trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
    name = strcat('caseDegenerated',num2str(i));
    title(name);
    exportgraphics(gca,strcat(name,'.png'),'Resolution',1000);

    hold off
end

keyboard;

figure
hold on
view(86,44);
axis ([-1,11,-1,11,-1,11]);
plotSolid(case1,conec + 1,'cyan','blue');
findInterior(A);
keyboard;
figure
hold on
view(86,44);
axis ([-1,11,-1,11,-1,11]);
case1_scale = transpose(0.95*eye(3) * transpose(case1 + [0.1,0.1,0.1]));
plotSolid(case1,conec + 1,'cyan','blue');
plotPlane(A)
%exportgraphics(gca,'halfplanes.png','Resolution',1000);
hold off


keyboard;
% for i = 1:5
%     pt_i = pts{i};
%     figure
%     hold on
%     view(30,30);    
%     plotSolid(pt_i,conec + ones(12,3),'cyan','blue');
%     title(strcat('caseDegenerated',num2str(i)));
%     hold off
% end
for i = 1:5
    pt_i = pts{i};
    
    pt_other = transform(pt_i);
    if (i == 2 || i == 3)
    pt_other = trans(pt_i);
    end
    figure
    hold on
    view(30,30);    
    plotSolid(pt_i,conec + ones(12,3),'cyan','blue');
    plotSolid(pt_other,conec + ones(12,3),'red','yellow');
    name = strcat('caseDegenerated',num2str(i));
    title(name);
    exportgraphics(gca,strcat(name,'.png'),'Resolution',1000);

    hold off
end

keyboard
pt_other = transform(pt,1);

figure
hold on
view(30,30);

plotSolid(pt,conec + ones(12,3),'cyan');
plotSolid(pt_other,conec + ones(12,3),'red');
hold off

caso2 = [ 0, 0, 8 ; 10, 0, 10; 0, 10, 10;10, 10, 8;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
pt_other = transform(caso2,1);


figure
hold on
view(30,30);

plotSolid(caso2,conec + ones(12,3),'cyan');
plotSolid(pt_other,conec + ones(12,3),'red');
hold off

fecho_pt_A = [42.8147392,	-94.2119827,-2.32802582	;...
44.5470543,	-93.2124481,-2.32802582	;...
41.8152008,	-92.4796677,-2.32802582	;...
43.5475159,	-91.4801254,-2.32802582	;...
42.8147392,	-94.2119827,-2.44060326	;...
44.5470543,	-93.2124481,-2.44060326	;...
41.8152008,	-92.4796677,-2.44060326	;...
43.5475159,	-91.4801254,-2.44060326	];


bbA = createBB(fecho_pt_A);
fecho_pt_B = [42.0000000,	-93.0000000	,-2.44060326;...	
44.0000000,	-93.0000000	,-2.44060326;...	
42.0000000,	-91.0000000	,-2.44060326;...	
44.0000000,	-91.0000000	,-2.44060326;...	
42.0000000,	-93.0000000	,-2.55318069;...	
44.0000000,	-93.0000000	,-2.55318069;...	
42.0000000,	-91.0000000	,-2.55318069;...	
44.0000000,	-91.0000000	,-2.55318069];	

media = sum(fecho_pt_B)/8;

[k1,av1] = convhull(fecho_pt_A);
n = [];
d = [];
figure
hold on
view(30,30)
axis equal
for i=1:size(k1,1)
    p0 = fecho_pt_A(k1(i,1),:);
    p1 = fecho_pt_A(k1(i,2),:);
    p2 = fecho_pt_A(k1(i,3),:);
    u = p1 - p0;
    v = p2 - p0;
    ni = cross(u,v)/ norm(cross(u,v)); 
    n(end+1,1:3) = ni;
    d(end+1) = dot(n(i,:),p0);
    %p = [p0;p1;p2];
    %trisurf([1,2,3],p(:,1),p(:,2),p(:,3),'FaceColor','cyan');
   % avg = (p0+p1+p2)/3;
%     quiver3(avg(1),avg(2),avg(3),100*ni(1),100*ni(2),100*ni(3));
%     quiver3(0,0,0,d(i)*ni(1),d(i)*ni(2),d(i)*ni(3));
%     drawPlan(ni',-d(i),'r');
end
[k2,av2] = convhull(fecho_pt_B);
for i=1:size(k2,1)
    p0 = fecho_pt_A(k2(i,1),:);
    p1 = fecho_pt_A(k2(i,2),:);
    p2 = fecho_pt_A(k2(i,3),:);
    u = p1 - p0;
    v = p2 - p0;
    ni = cross(u,v)/ norm(cross(u,v)); 
    n(end+1,1:3) = ni;
    d(end+1) = dot(n(i,:),p0);
    p = [p0;p1;p2];
    %trisurf([1,2,3],p(:,1),p(:,2),p(:,3),'FaceColor','cyan');
    %avg = (p0+p1+p2)/3;
     %quiver3(avg(1),avg(2),avg(3),100*ni(1),100*ni(2),100*ni(3));
     %quiver3(0,0,0,d(i)*ni(1),d(i)*ni(2),d(i)*ni(3));
     %drawPlan(ni',-d(i),'r');
end
hold off
center = [1159.28308,1240.57654,-3060.49463];
for i = 1:size(n,1)
    ni = n(i,1:3);
    di = d(i);    
    pos = ni*di;
    u = center - pos;
    %plot3(pos(1),pos(2),pos(3),'o');
    v = n;
    %drawPlan(n',d,'r');
    %quiver3(pos(1),pos(2),pos(3),1000*n(1),1000*n(2),1000*n(3));
    x = u / norm(u);
   % quiver3(pos(1),pos(2),pos(3),u(1),u(2),u(3));
    value = dot(x,ni);
    if (value > 0)
        break;
    end
end
figure
hold on
view(30,30);
h1 = trisurf(k1,fecho_pt_A(:,1),fecho_pt_A(:,2),fecho_pt_A(:,3),'FaceColor','cyan','FaceAlpha',0.5);
h2 = trisurf(k2,fecho_pt_B(:,1),fecho_pt_B(:,2),fecho_pt_B(:,3),'FaceColor','red','FaceAlpha',0.5);
% [kBbA,av2] = convhull(bbA);
%[kBbB,av2] = convhull(bbB);
%h1 = trisurf(kBbA,bbA(:,1),bbA(:,2),bbA(:,3),'FaceColor','g');
%h2 = trisurf(kBbB,bbB(:,1),bbB(:,2),bbB(:,3),'FaceColor','y');
axis equal 
A = [0.00000000   , 0.00000000  ,1.00000000  ;...
-0.00000000  , -0.00000000 ,-1.00000000 ;...
0.866159260  , 0.499768049 ,0.00000000  ;...
-0.866158307 , -0.499769688,-0.00000000 ;...
0.499768257  , -0.866159141,0.00000000  ;...
-0.499771088 , 0.866157472 ,0.00000000  ;...
0.00000000   , 0.00000000  ,1.00000000  ;...
0.00000000   , 0.00000000  ,-1.00000000 ;...
1.00000000   , 0.00000000  ,0.00000000  ;...
-1.00000000  , -0.00000000 ,-0.00000000 ;...
0.00000000   , -1.00000000 ,-0.00000000 ;...
0.00000000   , 1.00000000  ,0.00000000  ];

A(:,4) = 1.0;

b = [2.32802582 ;...
-2.44060326;...
7.99975967;...
-9.99995041;...
-103.000015;...
100.999985;...
2.44060326;...
-2.55318069;...
-44.0000000;...
42.0000000;...
-93.0000000;...
91.0000000];

% [A,b] = getAb(solid1,solid2);
%A = A(1:6,:);
%b = b(1:6,:);
n = size(A,2);
m = size(A,1);
c = zeros(n,1);
c(end) = -1;
center = linprog(c,A,-b);
r = center(4,1);
center = center(1:3);
%plot3(center(1),center(2),center(3),'o','MarkerSize',5,'MarkerFaceColor','r');
%plot3(1159.28308,1240.57654,-3060.49463,'o','MarkerSize',10,'MarkerFaceColor','g');
deslocated_point = center + (media' - center)*0.1;
plot3(43.547516406778271, -91.480126443935518, -2.4406032562255859,'o','MarkerSize',10,'MarkerFaceColor','g');

%plot3(1123.00610 ,1227.43628 ,-3050.73560,'o','MarkerSize',10,'MarkerFaceColor','b');
hold off

dual = [0.0000000000000000, 0.0000000000000000, 8.8827745116310588 ;...
-0.0000000000000000, -0.0000000000000000, -0.0000000000000000 ;...
0.0000000000000000, 0.0000000000000000, 0.0000000000000000 ;...
-0.43307931416594420, -0.24988493664565051, -0.0000000000000000 ;...
0.24988398660240710, -0.43307932448378439, 0.0000000000000000 ;...
0.0000000000000000, 0.0000000000000000, 0.0000000000000000 ;...
0.0000000000000000, 0.0000000000000000, -8.8827745116310588 ;...
2.2100248826259081, 0.0000000000000000, 0.0000000000000000 ;...
-0.64619670306557231, -0.0000000000000000, -0.0000000000000000 ;...
-0.0000000000000000, -0.65794946955283040, -0.0000000000000000 ;...
0.0000000000000000, 2.0827846760598385, 0.0000000000000000 ];
[kDual,av1] = convhull(dual);
figure
hold on
plot3(dual(:,1),dual(:,2),dual(:,3),'o');
view(30,30)
%trisurf(kDual,dual(:,1),dual(:,2),dual(:,3),'FaceColor','cyan');
hold off

finalPoints = zeros(size(kDual,1),3);
for i = 1:size(kDual,1)
    inc = kDual(i,:);
    p1 = dual(inc(1),:);
    p2 = dual(inc(2),:);
    p3 = dual(inc(3),:);
    normal = cross(p2 - p1, p3 - p1);
    normN  = norm(normal);
    normal = normal / norm(normal);
    d = dot(normal,p1);
    finalPoints(i,:) = (1/d)*normal + center';
end

[kDual,av1] = convhull(finalPoints);
figure
hold on
view(30,30)
trisurf(kDual,finalPoints(:,1),finalPoints(:,2),finalPoints(:,3),'FaceColor','green');
h1 = trisurf(k1,fecho_pt_A(:,1),fecho_pt_A(:,2),fecho_pt_A(:,3),'FaceColor','cyan','FaceAlpha',0.5);
h2 = trisurf(k2,fecho_pt_B(:,1),fecho_pt_B(:,2),fecho_pt_B(:,3),'FaceColor','red','FaceAlpha',0.5);
axis equal 

hold off

function drawPlan(n,d,color,A,xo)
if nargin == 1
    A = 20000;
    xo = [0;0;0];
    color = [1,0,0];
    d = 0;
elseif nargin == 2
    xo = [0;0;0];
    color = [1,0,0];
    A = 20000;
elseif nargin == 3
    xo = [0;0;0];
    A = 20000;
end
[x,y] = findTriedro(n);
xp1 = xo + A*x - n*d ;
yp1 = xo + A*y - n*d ;
xp2 = xo - A*x - n*d ;
yp2 = xo - A*y - n*d ;
h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],color);
set(h1, 'facealpha',0.6);
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
function bb = createBB(points)
bb = zeros(8,3);
pMax = [max(points(:,1)), max(points(:,2)), max(points(:,3))];
pMin = [min(points(:,1)), min(points(:,2)), min(points(:,3))];
bb(1,:) = pMin;
bb(2,:) = [pMin(1),pMax(2),pMin(3)];
bb(3,:) = [pMax(1),pMax(2),pMin(3)];
bb(4,:) = [pMax(1),pMin(2),pMin(3)];

bb(5,:) = pMax;
bb(6,:) = [pMin(1),pMax(2),pMax(3)];
bb(7,:) = [pMin(1),pMin(2),pMax(3)];
bb(8,:) = [pMax(1),pMin(2),pMax(3)];

end
function plotSolid(pt,conec,color,color_point)

for i = 1:size(pt,1)
pt_i = pt(i,:);
plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',10,'MarkerFaceColor',color_point);
text(pt_i(1),pt_i(2),pt_i(3),num2str(i), 'VerticalAlignment','top','HorizontalAlignment','left');
end
trisurf(conec,pt(:,1),pt(:,2),pt(:,3),'FaceAlpha',0.5,'FaceColor',color);
% for i = 1:12
%     id = conec(i,:);
%     trisurf(id,pt(:,1),pt(:,2),pt(:,3),'FaceAlpha',0.5,'FaceColor',color);
% end
end
function pt_final = transform(pt)
teta = 1;
centroide = sum(pt)/size(pt,1);
T = [cos(teta),-sin(teta),0;sin(teta),cos(teta),0;0,0,1];
pt_trans = pt - centroide;
pt_rot = transpose(T * pt_trans');
pt_final = pt_rot + centroide;
id = 3;
r = pt(id,:) - pt_final(id,:);
pt_final = pt_final + r;
end
function pt_final = trans(pt)
id_i = 1;
id_j = 3;
trans = pt(id_i,:) - pt(id_j,:);
pt_final = pt + trans;
end
function planes = BuildHalfPlanes(conec,pts)
planes = zeros(6,4);
for i = 1:2:11
    id = conec(i,:);
    pts_i = pts(id,:);
    [n_i,d_i] = getPlane(pts_i);
    id = conec(i+1,:);
    pts_i = pts(id,:);
    [n_j,d_j] = getPlane(pts_i);
    norm_ni = norm(n_i);
    norm_nj = norm(n_j);
    if (norm(n_i) == 0 || isnan(norm_ni))
        if (norm(n_j) == 0 || isnan(norm_nj))
            continue;
        else
            planes((i+1)/2,:) = [n_j,d_j];
        end
    else
        if (norm(n_j) == 0 || isnan(norm_nj))
            planes((i+1)/2,:) = [n_i,d_i];
        else
            n_medio = 0.5 * (n_i + n_j);
            d_medio = 0.5 * (d_i + d_j);
            planes((i+1)/2,:) = [n_medio,d_medio];
        end
    end
end

end
function [n,d] = getPlane(pts_i)
    u = pts_i(1,:) - pts_i(2,:);
    v = pts_i(1,:) - pts_i(3,:);
    n = cross(u,v)/norm(cross(u,v));
    d = dot(n,pts_i(1,:));
end
function out = getVectorizedPoints(conec,pts_i)
out = [];
pts_i_trans = pts_i';
for i = 1:6
    id = conec(i,:) + 1;   
    v = pts_i_trans(:,id);
    out(end+1,1:12) = v(:);
end
end
function plotPlane(A)
for i = 1:6
    n = A(i,1:3);
    d = A(i,4);
    drawPlan(n',d,'r');
end
end

function [k2,dualPoints,vol] = findInterior(M,trans)
%A = A(1:6,:);
%b = b(1:6,:);
A = M(:,1:3);
A(:,4) = 1.0;
b = -M(:,4);
n = size(A,2);
m = size(A,1);
c = zeros(n,1);
c(end) = -1;
tic
center = linprog(c,A,b);
r = center(4,1);
center = center(1:3)';
% center = [1,2,3];
coord_dual = [];
for i = 1:m
    n = A(i,1:3);
    d = b(i);
    new_d = d - dot(center,n);
    if (new_d == 0)
        continue;
    end
    coord_dual(end+1,1:3) = ( 1/new_d ) * n;
end
[k1,av1] = convhull(coord_dual);
dualPoints = zeros(size(k1,1),3);
count = 1;
for i = 1:size(k1,1)
    inc = k1(i,:);
    p1 = coord_dual(inc(1),:);
    p2 = coord_dual(inc(2),:);
    p3 = coord_dual(inc(3),:);
    normal = cross(p2 - p1, p3 - p1);
    normal = normal / norm(normal);
    d = dot(normal,p1);
    if (d == 0)
        continue;
    end
    dualPoints(count,:) = (1/d)* normal + center + trans;
    count = count + 1;
end
[k2,vol] = convhull(dualPoints);
% figure
% view(30,30);
% hold on
%trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
%hold off;
end
function volume = calculateVolume(cell,conec,id)
n = size(conec,1);
triangles = [];
for i = 1:n
    ids = conec(i,:);
    pts = cell(ids,:);
    u = pts(1,:) - pts(2,:);
    v = pts(1,:) - pts(3,:);
    area = 0.5*norm(cross(u,v));
    if area > 0
        triangles(end+1,1:3) = ids;
    end
end
if id == 4
    %could be concave
    volume = 0;    
    tetra_1 = triangles([2,3,4,7],:);
    volume = volume + getVol(cell,tetra_1);
    tetra_2 = triangles([1,5,6,8],:);
    volume = volume + getVol(cell,tetra_2);
    return;
end

volume = getVol(cell,triangles);

end
function volume = getVol(cell,triangles)
n_tri = size(triangles,1);
centroide = zeros(1,3);
for i = 1:n_tri
    ids = triangles(i,1:3);
    pts = cell(ids,:);
    c = (pts(1,:) + pts(2,:) + pts(3,:))/3;
    centroide = centroide + c;
end
centroide = centroide / n_tri;
% calculate volume
volume = 0;
for i = 1:n_tri
    ids = triangles(i,1:3);
    pts = cell(ids,:);
    u = pts(1,:) - pts(2,:);
    v = pts(1,:) - pts(3,:);
    n = cross(u,v);
    norm_n = norm(n);
    n = n / norm_n;
    area = norm_n;
    d = centroide -  pts(1,:);
    h = abs(dot(n,d));
    volume = volume + h*area*(1/6);
    showTetra([pts;centroide]);
end
end
function showTetra(pts)
tri1 = [1,2,3];
tri2 = [1,2,4];
tri3 = [2,3,4];
tri4 = [1,3,4];
plot3(pts(4,1),pts(4,2),pts(4,3),'o','MarkerSize', 20, 'MarkerFaceColor', 'red');
triangles = [tri1;tri2;tri3;tri4];
trisurf(triangles,pts(:,1),pts(:,2),pts(:,3),'FaceAlpha',0.5);
end
