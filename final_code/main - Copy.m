clc;
clear;
close all;

figure 
hold on
axis off
view(-70,43);
% dupla sobpreposta


pts_7912 = [-344.500000 2185.00000 -3147.56006 ; ...
-432.500000 2225.00000 -3158.59009 ; ...
-392.500000 2095.00000 -3138.87012 ; ...
-480.500000 2145.00000 -3156.04004 ; ...
-344.500000 2185.00000 -3147.56006 ; ...
-432.500000 2225.00000 -3158.59009 ; ...
-392.500000 2095.00000 -3145.87012 ; ...
-480.500000 2145.00000 -3163.04004 ];

pts_1122 = [675.500000 715.000000 -2983.62988 ; ...
587.500000 755.000000 -2987.16992 ; ...
627.500000 625.000000 -2983.01001 ; ...
539.500000 675.000000 -2988.05005 ; ...
675.500000 715.000000 -2990.62988 ; ...
587.500000 755.000000 -2994.16992 ; ...
627.500000 625.000000 -2990.01001 ; ...
539.500000 675.000000 -2995.05005 ];


pts_11647 = [-344.500000 2185.00000 -3147.56006 ; ...
-432.500000 2225.00000 -3158.59009 ; ...
-392.500000 2095.00000 -3145.87012 ; ...
-480.500000 2145.00000 -3163.04004 ; ...
-344.500000 2185.00000 -3147.56006 ; ...
-432.500000 2225.00000 -3158.59009 ; ...
-392.500000 2095.00000 -3152.87012 ; ...
-480.500000 2145.00000 -3166.83008 ];



pts_228 = [2870.50000 765.000000 -3221.79004 ; ...
2782.50000 815.000000 -3217.43994 ; ...
2822.50000 675.000000 -3213.40991 ; ...
2734.50000 725.000000 -3204.87012 ; ...
2870.50000 765.000000 -3228.79004 ; ...
2782.50000 815.000000 -3224.43994 ; ...
2822.50000 675.000000 -3220.40991 ; ...
2734.50000 725.000000 -3211.87012 ];

pts_4941 = [715.500000 575.000000 -2988.51001 ; ...
627.500000 625.000000 -2990.01001 ; ...
667.500000 485.000000 -2985.94995 ; ...
579.500000 535.000000 -2987.32007 ; ...
715.500000 575.000000 -2995.51001 ; ...
627.500000 625.000000 -2997.01001 ; ...
667.500000 485.000000 -2992.94995 ; ...
579.500000 535.000000 -2994.32007 ];

pts_3963 = [2870.50000 765.000000 -3228.79004 ; ...
2782.50000 815.000000 -3224.43994 ; ...
2822.50000 675.000000 -3220.40991 ; ...
2734.50000 725.000000 -3211.87012 ; ...
2870.50000 765.000000 -3235.79004 ; ...
2782.50000 815.000000 -3231.43994 ; ...
2822.50000 675.000000 -3227.40991 ; ...
2734.50000 725.000000 -3218.87012 ];

pts_1946 = [-663.500000  2015.00000  -3175.31006 ; ...
-751.500000  2065.00000  -3183.62988 ; ...
-711.500000  1925.00000  -3162.07007 ; ...
-799.500000  1975.00000  -3167.02002 ; ...
-663.500000  2015.00000  -3175.31006 ; ...
-751.500000  2065.00000  -3183.62988 ; ...
-711.500000  1925.00000  -3167.51001 ; ...
-799.500000  1975.00000  -3173.58008 ];

pts_15546 = [-615.500000   2105.00000  -3177.80005 ; ...
-703.500000   2145.00000  -3187.13989 ; ...
-663.500000   2015.00000  -3169.28003 ; ...
-751.500000   2065.00000  -3176.65991 ; ...
-615.500000   2105.00000  -3183.31006 ; ...
-703.500000   2145.00000  -3194.13989 ; ...
-663.500000   2015.00000  -3175.31006 ; ...
-751.500000   2065.00000  -3183.62988 ];

pts_11646 = [-432.500000 2225.00000 -3158.59009 ; ... 
-519.500000 2275.00000 -3170.91992 ; ... 
-480.500000 2145.00000 -3163.04004 ; ... 
-567.500000 2185.00000 -3178.47998 ; ... 
-432.500000 2225.00000 -3158.59009 ; ... 
-519.500000 2275.00000 -3170.91992 ; ... 
-480.500000 2145.00000 -3166.83008 ; ... 
-567.500000 2185.00000 -3180.59009 ];
pts_564 = [3029.50000 225.000000 -3199.95996 ; ...
2941.50000 265.000000 -3198.75000 ; ...
2981.50000 135.000000 -3192.76001 ; ...
2893.50000 185.000000 -3192.54004 ; ...
3029.50000 225.000000 -3206.95996 ; ...
2941.50000 265.000000 -3205.75000 ; ...
2981.50000 135.000000 -3199.76001 ; ...
2893.50000 185.000000 -3199.54004 ];

pts_4299 = [3029.50000 225.000000 -3206.95996 ; ...
2941.50000 265.000000 -3205.75000 ; ...
2981.50000 135.000000 -3199.76001 ; ...
2893.50000 185.000000 -3199.54004 ; ...
3029.50000 225.000000 -3213.95996 ; ...
2941.50000 265.000000 -3212.75000 ; ...
2981.50000 135.000000 -3206.76001 ; ...
2893.50000 185.000000 -3206.54004 ];

pts_572 = [3151.57812 -2008.25000 -3190.63989 ; ...
3100.35938 -1921.25000 -3192.01001 ; ...
3064.98438 -2062.25000 -3172.59009 ; ...
3013.79688 -1975.25000 -3175.23999 ; ...
3151.57812 -2008.25000 -3197.63989 ; ...
3100.35938 -1921.25000 -3199.01001 ; ...
3064.98438 -2062.25000 -3179.59009 ; ...
3013.79688 -1975.25000 -3182.23999 ];

pts_53067 = [3085.82812  -1954.75000  -3179.82007 ; ...
2997.82812  -1904.75000  -3174.13989 ; ...
3037.82812  -2044.75000  -3179.85010 ; ...
2949.82812  -1994.75000  -3173.11011 ; ...
3085.82812  -1954.75000  -3179.82007 ; ...
2997.82812  -1904.75000  -3174.13989 ; ...
3037.82812  -2044.75000  -3181.20996 ; ...
2949.82812  -1994.75000  -3174.10010 ];

pts_66031 = [-353.921875 -862.250000 -3088.35010 ; ...
-404.265625 -775.250000 -3090.58008 ; ...
-440.484375 -916.250000 -3090.69995 ; ...
-490.828125 -829.250000 -3088.29004 ; ...
-353.921875 -862.250000 -3088.35010 ; ...
-404.265625 -775.250000 -3092.51001 ; ...
-440.484375 -916.250000 -3090.69995 ; ...
-490.828125 -829.250000 -3088.29004 ];
pts_69093 = [-320.171875 -894.750000 -3079.85010 ; ...
-408.171875 -844.750000 -3084.42993 ; ...
-368.171875 -974.750000 -3087.70996 ; ...
-456.171875 -934.750000 -3092.04004 ; ...
-320.171875 -894.750000 -3086.85010 ; ...
-408.171875 -844.750000 -3091.42993 ; ...
-368.171875 -974.750000 -3094.70996 ; ...
-456.171875 -934.750000 -3099.04004 ];
pts_146 = [3007.82812  -805.250000  -3230.40991 ; ...
2957.48438  -718.750000  -3223.41992 ; ...
2921.26562  -859.250000  -3229.18994 ; ...
2870.92188  -772.750000  -3221.79004 ; ...
3007.82812  -805.250000  -3237.40991 ; ...
2957.48438  -718.750000  -3230.41992 ; ...
2921.26562  -859.250000  -3236.18994 ; ...
2870.92188  -772.750000  -3228.79004 ];

pts_397 = [2999.82812  -654.750000  -3227.76001 ; ...
2911.82812  -604.750000  -3219.10010 ; ...
2951.82812  -744.750000  -3213.13989 ; ...
2863.82812  -694.750000  -3206.43994 ; ...
2999.82812  -654.750000  -3234.76001 ; ...
2911.82812  -604.750000  -3226.10010 ; ...
2951.82812  -744.750000  -3220.13989 ; ...
2863.82812  -694.750000  -3213.43994 ];
pts_4132 = [2999.82812 -654.750000 -3234.76001 ; ...
2911.82812 -604.750000 -3226.10010 ; ...
2951.82812 -744.750000 -3220.13989 ; ...
2863.82812 -694.750000 -3213.43994 ; ...
2999.82812 -654.750000 -3241.76001 ; ...
2911.82812 -604.750000 -3233.10010 ; ...
2951.82812 -744.750000 -3227.13989 ; ...
2863.82812 -694.750000 -3220.43994 ];

pts_1 = pts_146;
pts_2 = pts_4132;

[cell_intersection_obj,volume] = CellIntersection(pts_1,pts_2,1);
keyboard;
bb = getBB([pts_1;pts_2]);
axis(bb);

[proj_coords_1, planes_1] = init(pts_1);
[proj_coords_2, planes_2] = init(pts_2);

figure 
hold on
axis off
view(-70,43);
axis(bb);

cube1 = Cube(proj_coords_1);
cube2 = Cube(proj_coords_2);
cube1.show('red');
cube2.show('blue');
[bool,p] = GJK(cube1,cube2);
plotPt(p);
bool_2 = check(p, [planes_1;planes_2]);
disp(bool);



function [proj_coords ,planes]  = init(pts)
showSurface(pts);
planes = CreateHalfs(pts);
plotPlaness(planes(1:2,:));
proj_coords = move(pts,planes);
end

function showSurface(pts)
ksi = linspace(0,1,10);
eta = linspace(0,1,10);
[ksi,eta] = meshgrid(ksi,eta);
planes = cell(6,1);
planes{1} = {pts(1,:),pts(2,:),pts(4,:), pts(3,:)};
planes{2} = {pts(5,:),pts(7,:),pts(8,:), pts(6,:)};
planes{3} = {pts(1,:),pts(3,:),pts(7,:), pts(5,:)};
planes{4} = {pts(2,:),pts(6,:),pts(8,:), pts(4,:)};
planes{5} = {pts(1,:),pts(5,:),pts(6,:), pts(2,:)};
planes{6} = {pts(3,:),pts(4,:),pts(8,:), pts(7,:)};
%% Draws
%axis off
%hold on
%view(-18,-23)
%text
for i = 1:size(pts,1)
    pt_i = pts(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','r');
    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i-1), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
end
k = [1,2,3;1,3,4;1,4,5;1,5,2];
R = sum(pts,1)/8;
plot3(R(1),R(2),R(3),'o','MarkerSize',10,'MarkerFaceColor','b');
for i = 1:6
    coords = planes{i};
    p0 = coords{1}';
    p1 = coords{2}';
    p2 = coords{3}';
    p3 = coords{4}';
    
    x = p0(1) + ksi.*(p1(1) - p0(1)) + eta.*(p3(1) - p0(1)) + ksi.*eta.*( (p2(1)-p0(1)) - (p1(1)-p0(1)) - (p3(1) - p0(1)) );
    y = p0(2) + ksi.*(p1(2) - p0(2)) + eta.*(p3(2) - p0(2)) + ksi.*eta.*((p2(2)-p0(2)) - (p1(2)-p0(2)) - (p3(2) - p0(2)));
    z = p0(3) + ksi.*(p1(3) - p0(3)) + eta.*(p3(3) - p0(3)) + ksi.*eta.*((p2(3)-p0(3)) - (p1(3)-p0(3)) - (p3(3) - p0(3)));
    surf(x,y,z,'FaceAlpha',0.5,'EdgeAlpha',0.1);
    pyramid = [R;p0';p1';p2';p3'];
    trisurf(k,pyramid(:,1),pyramid(:,2),pyramid(:,3),'FaceColor','green','FaceAlpha',0.2);
end
end
function planes = CreateHalfs(pts)
faces = cell(6,1);
faces(1) = {[ pts(1,:) ; pts(2,:) ; pts(4,:) ; pts(3,:)] };
faces(2) = {[ pts(5,:) ; pts(7,:) ; pts(8,:) ; pts(6,:)] };
faces(3) = {[ pts(2,:) ; pts(6,:) ; pts(8,:) ; pts(4,:)] };
faces(4) = {[ pts(1,:) ; pts(3,:) ; pts(7,:) ; pts(5,:)] };
faces(5) = {[ pts(1,:) ; pts(5,:) ; pts(6,:) ; pts(2,:)] };
faces(6) = {[ pts(3,:) ; pts(4,:) ; pts(7,:) ; pts(7,:)] };
planes = zeros(6,4);
for i = 1:6
    face_pts = faces{i};
    planes(i,:) = CreateFittingPlane(face_pts);
end
opp_ids = [2,1,4,3,6,5];
for i = 1:6
    plane = planes(i,:);
    normal = plane(1:3);
    %correct plane
    if (norm(normal) == 0)
        opp_plane = planes(opp_ids(i),:);
        normal = opp_plane(1:3);
        if (norm(normal) == 0)
            % não tem o que fazer
        else
            max_dist = - realmax;
            id_far = -1;
            plane_pts = faces{i};
            pn = - opp_plane(1:3) * opp_plane(4);
            for j = 1:4
                pi = plane_pts(j,:) - pn;
                dist = abs(dot(pi, normal));
                if dist > max_dist
                    id_far = j;
                    max_dist = dist;
                end
            end
            planes(i,:) = [-normal, dot(plane_pts(id_far,:),normal)];
        end
    end
end
end

function new_coords = move(pts,planes)
new_coords = zeros(8,3);

ids_pts = [1,2,4,3,5,6,8,7];
pts = pts(ids_pts,:);   % put on circular order

ids_planes = [4,5,3,6,1,2]; % lef front right back top bottom
planes = planes(ids_planes,:);
view(62,80);

[o,d] = intersect3D(-planes(6,4)*planes(6,1:3),-planes(5,4)*planes(5,1:3));

for i = 1:4
    p_i = pts(i,:);
    p_j = pts(i+4,:);
    edge = p_i - p_j;
    if (dot(edge,edge) == 0)
        edge_dir = cross(planes(i,:),planes(mod(i,4)+1,:));
        p_i = p_j + edge_dir;
        p_j = p_j - edge_dir;
        edge = p_i - p_j;
    end
    p_i = rayPlaneIntersection(p_i,edge,planes(5,:));
    p_j = rayPlaneIntersection(p_j,edge,planes(6,:));
    new_edge = p_i - p_j;
    if dot(edge,new_edge) < 0
        p_j = rayPlaneIntersection(o,d,planes(i,:));
        p_i = rayPlaneIntersection(o,d,planes(mod(i,4)+1,:));
    end
    new_coords(i,:) = p_i ;
    new_coords(i + 4,:) = p_j ;
end
end
function [o,z] = intersect3D(n,m)
z = cross(m,n);
m_perp = cross(z,m);
c = dot(n,n);
d = dot(m,n);
e = dot(m_perp,n);
lam2 = (c-d)/e;
o = m + lam2*m_perp;
lam = 10000;
p1 = o + lam*z/norm(z);
p2 = o - lam*z/norm(z);
d = z/norm(z);
line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'linewidth',3,'color','black');
check1 = dot(o-n,n);
check2 = dot(o-m,m);
end
function p = rayPlaneIntersection(o,ray,n)
d = n(4);
n_vec = n(1:3);
d_i = dot(n_vec,o) + d;
lam = - d_i / dot(n_vec,ray);
p = o + lam*ray;
end
function plane = CreateFittingPlane(pts)
centroide = mean(pts);
normal = zeros(1,3);
plane = zeros(1,4);
%k = convhull(pts(:,1),pts(:,2),pts(:,3));

%a = getNormal(pts([1,2,3],:));
% b = getNormal(pts([1,3,4],:));
% if dot(a,b) < 0
%     b = -b;
% end
% normal = a + b;
% for i = 1:size(k,2)
%     normal = normal + getNormal(pts(k(i,:),:));
% end
for i = 1:4
    p0 = pts(i,:);
    p1 = pts(mod(i,4) + 1 , :);
    p2 = pts(mod(i+1,4) + 1 , :);
    normal_i = cross(p1-p0,p2-p0);
    len = norm(normal_i);
    if len == 0
        continue
    end
    normal_i = normal_i / len;
    normal = normal + normal_i;
end
if norm(normal) == 0
    return
end
normal = normal/norm(normal);
plane = [normal,-dot(normal, centroide)];
end
function bb = getBB(coords)
margin = 10.0;
x_min = min(coords(:,1)) - margin;
x_max = max(coords(:,1)) + margin;
y_min = min(coords(:,2)) - margin;
y_max = max(coords(:,2)) + margin;
z_min = min(coords(:,3)) - margin;
z_max = max(coords(:,3)) + margin;
bb = [x_min,x_max,y_min,y_max,z_min,z_max];
end
function plotPlaness(planes)
for i = 1:size(planes,1)
    n = planes(i,1:3);
    d = planes(i,4);
    drawPlan(n',d,'r');
end
end
function h1 = drawPlan(n,d,color,A,xf,xo)
if nargin == 1
    A = 20000;
    xo = [0;0;0];
    color = [1,0,0];
    d = 0;
    xf = -n*d;
elseif nargin == 2
    xo = [0;0;0];
    color = [1,0,0];
    A = 20000;
    xf = -n*d;
elseif nargin == 3
    xo = [0;0;0];
    A = 2000000;
    xf = -n*d;
elseif nargin == 5
    xo = [0;0;0];
end
[x,y] = findTriedro(n);
xp1 = xo + A*x + xf ;
yp1 = xo + A*0.1*y + xf ;
xp2 = xo - A*x + xf;
yp2 = xo - A*0.1*y + xf;
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
function [bool,p] = GJK(convex1,convex2)

plotPt([0,0,0],'yellow');

triA = Tetrahedra();
triB = Tetrahedra();
tri = Tetrahedra();
d0 = (convex2.centroide - convex1.centroide);
d0 = d0/norm(d0);
d0 = [0.4,0.7,0.6];
[vA,vB] = support(convex1,convex2,d0);
A = vA - vB;
%plot(vA(1),vA(2),'*');
%plot(vB(1),vB(2),'*');
if dot(A,A) == 0
    p = vA;
    bool = true;
    return
end
plotPt(A);
tri.append(A);
triA.append(vA);
triB.append(vB);
d = -A;
d = d/norm(d);
while true
    [vA,vB] = support(convex1,convex2,d);
    P = vA - vB;
    plotPt(P,'b');
    dot_prod = (dot(P,d)/(norm(P)*norm(d)))
    dotp = dot(P,d)
    if dot(P,d) <= 1e-7
        bool = false;
        break
    end
    tri.append(P);
    triA.append(vA);
    triB.append(vB);
    [flag,d] = handleSimplex(tri,triA,triB);
    if flag
        bool = true;
        break;
    end
end
lam = tri.getOriginBary();
% figure
% view(30,30)
% hold on
% triA.plotTetra();
% convex1.show();
%triB.plotTetra();
dir_A = triA.checkDegenerate();
dir_B = triB.checkDegenerate();
p = triA.plotBaryPoint(lam);
%plotPt(p);
%p = triA.pertubate(p,convex1.pts,[0.0744385421 , 0.0297401678 , 0.996782064 ],lam);
%p = triB.pertubate(p,convex1.pts,[0.0744385421 , 0.0297401678 , 0.996782064 ],lam);
%p = triB.plotBaryPoint(lam);
%p = p + 0.5*dir_B;

end

function plotPt(P,color)
if nargin == 1 
    color = 'black';
end
plot3(P(1),P(2),P(3),'o','MarkerFaceColor',color,'MarkerSize',5);
end
function [vA,vB] = support(convex1,convex2,di)
vA = convex1.suportfunction(di');
vB = convex2.suportfunction(-di');
v =  vA - vB;
d_unit = 0.4*di/norm(di);
% quiver3(convex1.centroide(1),convex1.centroide(2),convex1.centroide(3),d_unit(1),d_unit(2),d_unit(3),'color','magenta');
% quiver3(convex2.centroide(1),convex2.centroide(2),convex2.centroide(3),-d_unit(1),-d_unit(2),-d_unit(3),'color','magenta');
% plot3(vA(1),vA(2),vA(3),'o');
% plot3(vB(1),vB(2),vB(3),'o');
end
function [bool, d] = handleSimplex(tetra,tetraA,tetraB)
if tetra.n == 2
    [bool,d] = tetra.GetPerpDir2Line();
    return
elseif tetra.n == 3
    [bool,d] = tetra.GetPerpDir2Tri();
    return
end
[bool, d] = tetra.ContainOrigin(tetraA,tetraB);
end
function out = check(center,planes)
out = true;
for i = 1:size(planes,1)
    n = planes(i,1:3);
    d = planes(i,4);
    u = center + n * d;
    u = u / norm(u) ;
    angle = dot(u,n);
    if (angle >= 0.0)
         out = false;
    end
end
end