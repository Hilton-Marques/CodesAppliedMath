% Clear workspace
clear
close(findall(0,'Type','figure'));
clc

%convex_study()

v1 = [-0.5,-0.5,-0.5];
v2 = [0.5,-0.5,-0.5];
v3 = [0.5,0.5,-0.5];
v4 = [-0.5,0.5,-0.5];

v5 = [-0.5,-0.5,0.5];
v6 = [0.5,-0.5,0.5];
v7 = [0.5,0.5,0.5];
v8 = [-0.5,0.5,0.5];

%plotPlanScene();
%showCell();
% 
% 
% keyboard;
cell_146 = [3007.82812 , -805.250000 , -3230.40991 ; ...
2957.48438 , -718.750000 , -3223.41992 ; ...
2921.26562 , -859.250000 , -3229.18994 ; ...
2870.92188 , -772.750000 , -3221.79004 ; ...
3007.82812 , -805.250000 , -3237.40991 ; ...
2957.48438 , -718.750000 , -3230.41992 ; ...
2921.26562 , -859.250000 , -3236.18994 ; ...
2870.92188 , -772.750000 , -3228.79004 ];

cell_397 = [2999.82812, -654.750000 , -3227.76001 ; ...
2911.82812, -604.750000 , -3219.10010 ; ...
2951.82812, -744.750000 , -3213.13989 ; ...
2863.82812, -694.750000 , -3206.43994 ; ...
2999.82812, -654.750000 , -3234.76001 ; ...
2911.82812, -604.750000 , -3226.10010 ; ...
2951.82812, -744.750000 , -3220.13989 ; ...
2863.82812, -694.750000 , -3213.43994 ];


cell_4132 = [2999.82812 , -654.750000 , -3234.76001 ; ...
2911.82812 , -604.750000 , -3226.10010 ; ...
2951.82812 , -744.750000 , -3220.13989 ; ...
2863.82812 , -694.750000 , -3213.43994 ; ...
2999.82812 , -654.750000 , -3241.76001 ; ...
2911.82812 , -604.750000 , -3233.10010 ; ...
2951.82812 , -744.750000 , -3227.13989 ; ...
2863.82812 , -694.750000 , -3220.43994 ];

cell_146_b = [2999.82812 , -654.750000 , -3227.76001 ;...
2911.82812 , -604.750000 , -3219.10010 ;...
2951.82812 , -744.750000 , -3213.13989 ;...
2863.82812 , -694.750000 , -3206.43994 ;...
2999.82812 , -654.750000 , -3234.76001 ;...
2911.82812 , -604.750000 , -3226.10010 ;...
2951.82812 , -744.750000 , -3220.13989 ;...
2863.82812 , -694.750000 , -3213.43994 ];

cell_1 = [2804.32812 , -1045.75000 , -3219.10010 ;...
2753.14062 , -958.750000 , -3208.70996 ;...
2717.76562 , -1099.75000 , -3206.43994 ;...
2666.54688 , -1012.75000 , -3199.65991 ;...
2804.32812 , -1045.75000 , -3226.10010 ;...
2753.14062 , -958.750000 , -3215.70996 ;...
2717.76562 , -1099.75000 , -3213.43994 ;...
2666.54688 , -1012.75000 , -3206.65991 ];

cell_2 = [2807.82812 , -1004.75000 , -3217.83008 ; ...
2719.82812 , -954.750000 , -3211.27002 ; ...
2759.82812 , -1094.75000 , -3210.07007 ; ...
2671.82812 , -1044.75000 , -3202.40991 ; ...
2807.82812 , -1004.75000 , -3224.83008 ; ...
2719.82812 , -954.750000 , -3218.27002 ; ...
2759.82812 , -1094.75000 , -3217.07007 ; ...
2671.82812 , -1044.75000 , -3209.40991 ];

cell_391 = [2545.01562 , -619.250000 , -3164.29004 ; ...
2493.79688 , -531.750000 , -3144.40991 ; ...
2458.45312 , -673.250000 , -3146.87988 ; ...
2407.23438 , -585.750000 , -3126.55005 ; ...
2545.01562 , -619.250000 , -3171.29004 ; ...
2493.79688 , -531.750000 , -3151.40991 ; ...
2458.45312 , -673.250000 , -3153.87988 ; ...
2407.23438 , -585.750000 , -3133.55005 ];

cell_559 = [2552.82812 , -634.750000 , -3167.52002 ;...
2464.82812 , -594.750000 , -3154.93994 ;...
2504.82812 , -724.750000 , -3157.64990 ;...
2416.82812 , -674.750000 , -3141.03003 ;...
2552.82812 , -634.750000 , -3174.52002 ;...
2464.82812 , -594.750000 , -3161.93994 ;...
2504.82812 , -724.750000 , -3164.64990 ;...
2416.82812 , -674.750000 , -3148.03003];

cell_4796 = [2172.98438 , -1547.25000 , -3156.75000 ;...
2121.79688 , -1460.25000 , -3142.84009 ;...
2086.42188 , -1601.25000 , -3140.16992 ;...
2035.20312 , -1514.25000 , -3139.38989 ;...
2172.98438 , -1547.25000 , -3163.75000 ;...
2121.79688 , -1460.25000 , -3149.84009 ;...
2086.42188 , -1601.25000 , -3147.16992 ;...
2035.20312 , -1514.25000 , -3146.38989 ];

cell_46127 = [2208.82812 , -1474.75000 , -3164.09009 ;...
2120.82812 , -1424.75000 , -3159.55005 ;...
2160.82812 , -1564.75000 , -3158.36011 ;...
2072.82812 , -1514.75000 , -3152.73999 ;...
2208.82812 , -1474.75000 , -3171.09009 ;...
2120.82812 , -1424.75000 , -3166.55005 ;...
2160.82812 , -1564.75000 , -3165.36011 ;...
2072.82812 , -1514.75000 , -3159.73999 ];

cell_314 = [2942.10938 , -1079.25000 , -3244.08008 ; ...
2890.89062 , -991.750000 , -3231.48999 ; ...
2855.54688 , -1133.25000 , -3227.76001 ; ...
2804.32812 , -1045.75000 , -3219.10010 ; ...
2942.10938 , -1079.25000 , -3251.08008 ; ...
2890.89062 , -991.750000 , -3238.48999 ; ...
2855.54688 , -1133.25000 , -3234.76001 ; ...
2804.32812 , -1045.75000 , -3226.10010 ];

cell_19405 = [2895.82812 -1054.75000 -3221.62012 ; ...
2807.82812 -1004.75000 -3217.83008 ; ...
2847.82812 -1144.75000 -3213.53003 ; ...
2759.82812 -1094.75000 -3210.07007 ; ...
2895.82812 -1054.75000 -3228.62012 ; ...
2807.82812 -1004.75000 -3224.83008 ; ...
2847.82812 -1144.75000 -3220.53003 ; ...
2759.82812 -1094.75000 -3217.07007];

cell_42785 = [191.234375 , 174.250000 , -3017.79004 ; ...
135.046875 , 252.750000 , -3018.12012 ; ...
104.671875 , 120.250000 , -3017.87012 ; ...
53.4531250 , 207.250000 , -3020.08008 ; ...
191.234375 , 174.250000 , -3024.79004 ; ...
135.046875 , 252.750000 , -3025.12012 ; ...
104.671875 , 120.250000 , -3024.87012 ; ...
53.4531250 , 207.250000 , -3027.08008 ];

cell_949 = [119.828125  , 125.250000  , -3025.01001; ...
31.8281250  , 175.250000  , -3029.66992; ...
71.8281250  , 35.2500000  , -3020.76001; ...
-16.1718750 , 85.2500000 , -3028.15991; ...
119.828125  , 125.250000  , -3032.01001; ...
31.8281250  , 175.250000  , -3036.66992; ...
71.8281250  , 35.2500000  , -3027.76001; ...
-16.1718750 , 85.2500000  , -3035.15991];

cell_934737 = [1800.00000 , 6100.00000 , -5790.39990 ; ...
1900.00000 , 6100.00000 , -5772.85986 ; ...
1800.00000 , 6200.00000 , -5803.02002 ; ...
1900.00000 , 6200.00000 , -5784.47021 ; ...
1800.00000 , 6100.00000 , -5796.93018 ; ...
1900.00000 , 6100.00000 , -5779.50000 ; ...
1800.00000 , 6200.00000 , -5809.58984 ; ...
1900.00000 , 6200.00000 , -5791.14990 ];

cell_11497458 = [1844.12744 , 6127.64062 , -5780.93018 ; ...
1799.32007 , 6217.04053 , -5791.16016 ; ...
1754.72778 , 6082.83350 , -5788.06982 ; ...
1709.92041 , 6172.23291 , -5797.75000 ; ...
1844.12744 , 6127.64062 , -5785.04004 ; ...
1799.32007 , 6217.04053 , -5794.97998 ; ...
1754.72778 , 6082.83350 , -5792.12012 ; ...
1709.92041 , 6172.23291 , -5801.49023 ];

cell_16012 = [2743.70312 , -2488.75000 , -3157.94995; ...
2692.48438 , -2401.75000 , -3160.77002; ...
2657.14062 , -2542.75000 , -3142.62988; ...
2605.92188 , -2455.75000 , -3148.57007; ...
2743.70312 , -2488.75000 , -3164.94995; ...
2692.48438 , -2401.75000 , -3167.77002; ...
2657.14062 , -2542.75000 , -3149.62988; ...
2605.92188 , -2455.75000 , -3155.57007];

cell_72779 = [2709.82812 , -2434.75000 , -3151.72998 ; ...
2622.82812 , -2384.75000 , -3152.83008 ; ...
2661.82812 , -2524.75000 , -3149.22998 ; ...
2574.82812 , -2474.75000 , -3148.01001 ; ...
2709.82812 , -2434.75000 , -3151.72998 ; ...
2622.82812 , -2384.75000 , -3159.83008 ; ...
2661.82812 , -2524.75000 , -3156.22998 ; ...
2574.82812 , -2474.75000 , -3155.01001 ];

cell_11603  = [1200.00000 , -14300.0000 , -6057.75000 ; ...
1300.00000 , -14300.0000 , -6057.18994 ; ...
1200.00000 , -14200.0000 , -6054.35010 ; ...
1300.00000 , -14200.0000 , -6054.47998 ; ...
1200.00000 , -14300.0000 , -6064.37012 ; ...
1300.00000 , -14300.0000 , -6063.81006 ; ...
1200.00000 , -14200.0000 , -6060.95020 ; ...
1300.00000 , -14200.0000 , -6061.08008 ];

cell_133048 = [1100.00000 , -14400.0000 , -6067.66992 ; ...
1200.00000 , -14400.0000 , -6066.93018 ; ...
1100.00000 , -14300.0000 , -6064.56006 ; ...
1200.00000 , -14300.0000 , -6064.37012 ; ...
1100.00000 , -14400.0000 , -6074.31982 ; ...
1200.00000 , -14400.0000 , -6073.58984 ; ...
1100.00000 , -14300.0000 , -6071.18018 ; ...
1200.00000 , -14300.0000 , -6071.00000 ];

cell_trans_1 = [3007.82812 , -805.250000 , -3230.40991 ; ...
2957.48438 , -718.750000 , -3223.41992 ; ...
2870.92188 , -772.750000 , -3221.79004 ; ...
2921.26562 , -859.250000 , -3229.18994 ; ...
3007.82812 , -805.250000 , -3237.40991 ; ...
2957.48438 , -718.750000 , -3230.41992 ; ...
2870.92188 , -772.750000 , -3228.79004 ; ...
2921.26562 , -859.250000 , -3236.18994 ];
planes_1 = [0.0500084609 , -0.0538479611,  0.997296154   , 3027.99976  ; ...  
-0.0500084609, 0.0538479611 ,  -0.997296154  , -3034.98047 ; ...
-0.529282928 , 0.848445475  , -0.00000000  , 2175.16626  ; ...
0.529282928  , -0.848445475 ,  0.00000000  , -2275.20288 ; ...
0.864276767  ,  0.503016591 ,  -0.00000 , -2194.54199 ; ...
-0.864276767 , -0.503016591 , 0.00000000 ,  2092.56494 ]; 

%cell_trans_1 = move(cell_trans_1, planes_1);

cell_trans_2 = [2999.82812 , -654.750000 , -3227.76001 ; ...
2911.82812 , -604.750000 , -3219.10010 ; ...
2863.82812 , -694.750000 , -3206.43994 ; ...
2951.82812 , -744.750000 , -3213.13989 ; ...
2999.82812 , -654.750000 , -3234.76001 ; ...
2911.82812 , -604.750000 , -3226.10010 ; ...
2863.82812 , -694.750000 , -3213.43994 ; ...
2951.82812 , -744.750000 , -3220.13989];

planes_2 = [0.131468609 , 0.0796258524  ,0.988117337     ,2846.67236  ;...
-0.131468624, -0.0796258673 ,-0.988117337    ,-2853.58887 ;...
-0.882352948, 0.470588237   ,-0.00000000     ,2853.84839  ;...
0.882352948 , -0.470588237  ,0.00000000      ,-2955.02490 ;... 
0.494009405 , 0.869456530   ,-0.00000000     ,-912.666687 ;...  
-0.494009405, -0.869456530  ,0.00000000      ,810.70306] ;

%cell_trans_2 = move(cell_trans_2, planes_2);

planes  = [planes_1;planes_2];

planes = [-0.0767133608 , 0.0976413265  , 0.992260635  ,  2899.79248 ; ...
0.112258554   , -0.0449573807 , -0.992661476 ,  -2995.77710  ; ...
-0.882352948  , 0.470588237   , -0.00000000  ,  -1630.14709  ; ...
0.857492983   , -0.514495790  , 0.00000000   ,  1610.80054   ; ... 
0.413802952   , 0.910366535   , -0.00000000  ,  -1754.31775  ; ...
-0.413802952  , -0.910366535  , 0.00000000   ,  1661.62585  ];

pts_1 = [-485.734375 , -1998.25000 , -3043.43994 ; ...
-536.953125 , -1911.25000 , -3041.27002 ; ...
-623.515625 , -1965.25000 , -3049.64990 ; ...
-567.296875 , -2043.75000 , -3044.73999 ; ...
-485.734375 , -1998.25000 , -3043.43994 ; ...
-536.953125 , -1911.25000 , -3041.27002 ; ...
-623.515625 , -1965.25000 , -3054.12988 ; ...
-567.296875 , -2043.75000 , -3051.73999 ];

planes_1 = [-0.0476760305 , -0.0134470640 , 0.998772442  , 2988.03442  ; ...
0.0904662684  , 0.0568772927  , -0.994274020 , -2867.53687 ; ...
-0.529282928  , 0.848445415   , -0.00000000  , 1337.39111  ; ...
0.487176329   , -0.873303711  , 0.00000000   , -1508.44080 ; ...
0.813011229   , 0.582248092   , -0.00000000  , 1558.38477  ; ...
-0.813011229 ,  -0.582248092 ,  0.00000000   , -1651.18823 ];

projected_pp1 = convex_study(pts_1,planes_1);


pts_2 = [-449.171875 , -1964.75000 , -3037.35010 ; ...
-537.171875 , -1914.75000 , -3031.93994 ; ...
-585.171875 , -2004.75000 , -3042.62988 ; ...
-497.171875 , -2044.75000 , -3044.41992 ; ...
-449.171875 , -1964.75000 , -3044.35010 ; ...
-537.171875 , -1914.75000 , -3038.93994 ; ...
-585.171875 , -2004.75000 , -3049.62988 ; ...
-497.171875 , -2044.75000 , -3044.41992 ];

planes_2 = [-0.00966911297 , -0.0977790058 , 0.995161176  , 2825.55640  ; ...
0.0241735037   , 0.0466740616  , -0.998617589 , -2935.10498 ; ...
-0.882352948   , 0.470588237   , -0.00000000  , 427.083618  ; ... 
0.857492924    , -0.514495790  , 0.00000000   , -625.693848 ; ...
0.494009405    , 0.869456530   , -0.00000000  , 1930.15991  ; ...
-0.413802952   , -0.910366535  , 0.00000000   , -2067.20312 ];

projected_pp2 = convex_study(pts_2,planes_2);

planes = [planes_1;planes_2];
cube1 = Cube(projected_pp1);
cube2 = Cube(projected_pp2);
%cube1 = Cube([v1;v2;v3;v4;v5;v6;v7;v8]);
%cube2 = Cube([v1;v2;v3;v4;v5;v6;v7;v8]);
%cube2.transform(pi/3.45,0.6*[1.5,0,1.5]);
%cube1.transform(0,0.6*[1.5,0,1.5]);

figure
hold on
lighting gouraud;
a1 = camlight('right');
a2 = camlight('left');
a3 = camlight('headlight');
axis equal;
set(gca,'visible','off');
set(gcf,'color','white');
view(-28,17);

cube1.plot_tri();
cube2.plot_tri();
%exportgraphics(gca,'foto1.jpeg','Resolution',1000);


p_m = MikSub(cube1,cube2);
k = convhull(p_m(:,1),p_m(:,2),p_m(:,3));
trisurf(k,p_m(:,1),p_m(:,2),p_m(:,3),'FaceColor','None','EdgeAlpha',0.3);
%plot3(p_m(k,1),p_m(k,2),p_m(k,3),'*');
[bool,p] = GJK(cube1,cube2);
out = check(p,planes);
hold off
figure
hold on
lighting gouraud;
a1 = camlight('right');
a2 = camlight('left');
a3 = camlight('headlight');
axis equal;
set(gca,'visible','off');
set(gcf,'color','white');
view(-28,17);
plotPt(p,'yellow');
cube1.plot_tri();
cube2.plot_tri();
%plotPt(p,'red');
exportgraphics(gca,'GJK.jpeg','Resolution',1000)
findInterior(p);

%p_cheby = [2714.5909296287164, -1065.0591676686797, -3209.7786412814116];
%plotPt(p_cheby,'blue');

hold off
%%exportgraphics(gca,'foto2.jpeg','Resolution',1000);

function p_m = MikSub(convex1,convex2)
n1 = size(convex1.pts,1);
n2 = size(convex2.pts,1);
p_m = zeros(n1*n2,3);
count = 1;
for i = 1:n1
    pi = convex1.pts(i,:);
    for j = 1:n2
        pj = convex2.pts(j,:);
        p_m(count,:) = pi - pj;
        count = count + 1;
    end
end
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
    if dot(P,d) < 0
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
figure
view(30,30)
hold on
triA.plotTetra();
convex1.plot_tri();
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
function plotPt(P,color)
if nargin == 1 
    color = 'black';
end
plot3(P(1),P(2),P(3),'o','MarkerFaceColor',color,'MarkerSize',10);
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
function [k2,dualPoints,vol,bool,center] = findInterior(gjk)
M = [0.0726542175  -0.0163173564  0.997223616 , 2892.33643 ; ...
-0.0726542249 0.0163173564   -0.997223616 , -2899.31714 ; ...
-0.529282928  0.848445475    0.00000000 , 1160.61292 ; ...
0.529282928   -0.848445475   0.00000000 , -1260.64941 ; ...
0.864276767   0.503016591    0.00000000 , -671.937988 ; ...
-0.864276767  -0.503016591   0.00000000 , 569.961121 ; ...
0.0816126913  0.0311308019   0.996177852 , 2919.93262 ; ...
-0.0816126913 -0.0311308019  -0.996177852 , -2926.90601 ; ...
-0.882352948  0.470588237    0.00000000 , 1244.14246 ; ...
0.882352948   -0.470588237   0.00000000 , -1344.43652 ; ...
0.498283893   0.867013931    0.00000000 , 5.10083008 ; ...
-0.498283893  -0.867013931   0.00000000 , -107.049683 ];


trans = [0,0,0];
flag = true;
bool = false;
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
plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','green');
if flag
    x = GetChebyshevCenter(M);
end
%center = [2744.4892256379862, -1033.7425099725353, -3215.2021284151992];
%center = gjk;
gjk = [1119.4609678565462, -701.05068103115627, -3001.0229556558356];
err = M(:,1:3) * gjk' + M(:,4);

gjk = [1119.4609678565462 ,-701.05068103115627,-3001.0229556558356];
dir = center - gjk;
err = M(:,1:3) * gjk' + M(:,4);
vec = -M(2,1:3);
parallel = vec - dot(M(7,1:3),vec)*M(7,1:3);
parallel = parallel / norm(parallel);
dots = M(:,1:3) * parallel';

gjk = gjk + 44*parallel;
err = M(:,1:3) * gjk' + M(:,4);

last = cross(M(7,1:3),M(10,1:3));
dots = M(:,1:3) * last';
gjk = gjk - 42*last;
err = M(:,1:3) * gjk' + M(:,4);

parallel = gjk - dot(M(10,1:3),gjk)*M(10,1:3);
parallel = parallel / norm(parallel);
err = M(:,1:3) * gjk' + M(:,4);
dots = M(:,1:3) * -parallel';
gjk = gjk + 6.8*parallel;


for i = 1:size(M,1)
    n = M(i,1:3);
    d = M(i,4);
    u = center + n * d;
    u = u / norm(u) ;
    angle = dot(u,n);
    if (angle >= 0.0)
        k2 = 0;
        dualPoints = [];
        vol = 0;
        return;
    end
    
end
bool = true;
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
coord_points = coord_dual(unique(k1(:),'stable'),:);
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
    dualPoints(count,:) = (1/d)* normal + center + trans ;
    count = count + 1;
end
[k2,vol] = convhull(dualPoints);
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3));
end

function center = GetChebyshevCenter(M)
A = M(:,1:3);
A(:,4) = 1.0;
b = -M(:,4);
n = size(A,2);
m = size(A,1);
c = zeros(n,1);
c(end) = -1;
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
%Part 1
[x,points,cost,basis] = solveLP(newAA,b,newff,basis,b);
%Part 2
init = x(basis);
newf = zeros(20,1);
n = size(newf,1);
newAA(:,end-nIds+1:end) = [];
newf(4,1) = -1.0;
newf(8,1) = 1.0;
[x,points,cost,basis] = solveLP(newAA,b,newf,basis,init);

center = x(1:3)' - x(5:7)';

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
m = size(A,1);
points(end+1,1:3) = x(1:3,1) - x(4:6,1);
for iter = 1:100
    d = zeros(size(c,1),1);
    %dAux = zeros(size(c,1),1);
    y = Binv' * c(basis);           % this y may not be feasible
    nonBasics = setdiff(1:n,basis);
    rmin = 0;
    in = -1;
    for k = nonBasics
        rmin_in = c(k) - At(k,:)*y;
        if (rmin_in < -.00000001 && rmin_in < rmin)
            %break
            rmin = rmin_in;
            in = k;
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
    A*dAux;
    xb = x(basis);
    db = d(basis);
    teta = realmax;
    out = -1;
    for i = 1:m
        dbi = db(i);
        if (dbi > 1e-3 )
            inv_dbi = 1 / dbi;
            xbi = xb(i);
            value = xbi * inv_dbi;
            if value < teta
                teta = value;
                out = i;
            end
        end
    end
    if (out == -1)
        break;
    end
    
    cost = cost + teta*rmin  % lower cost at end of step
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
function plotPlanScene()
figure
hold on
lighting gouraud;
a1 = camlight('right');
a2 = camlight('left');
a3 = camlight('headlight');
axis equal;
set(gca,'visible','off');
set(gcf,'color','white');
view(-136,4);
n1 = [1;1;1];
n1 = n1 / norm(n1);
drawPlan(n1,10,'r');
n2 = [0;0;1];
drawPlan(n2,30,'b');
p = [2000;0;-5000];
plotPt(p);
%exportgraphics(gca,'planScene1.jpeg','Resolution',1000);
quiver3(0,0,0,0,0,10000);
%exportgraphics(gca,'planScene2.jpeg','Resolution',1000);
n3 = n2 - dot(n2,n1)*n1;
n3 = 10000*n3;
quiver3(0,0,0,n3(1),n3(2),n3(3),'color','black');
%exportgraphics(gca,'planScene3.jpeg','Resolution',1000);
quiver3(p(1),p(2),p(3),n3(1),n3(2),n3(3),'color','black');
%exportgraphics(gca,'planScene4.jpeg','Resolution',1000);
k = p + n3;
plotPt(k,'b');
exportgraphics(gca,'planScene5.jpeg','Resolution',1000);

end
function  h1 = drawPlan(n,d,color,A,xo)
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
function showCell()
ksi = linspace(0,1,10);
eta = linspace(0,1,10);
[ksi,eta] = meshgrid(ksi,eta);
pts = [0, 0, 10;...
    10, 0, 10;...
    0, 10, 10;...
    10, 10, 10;...
    0, 0, 0;...
    10, 0, 0;...
    0, 10, 0;...
    10, 10, 0];
lam1 = 1.8;
lam2 = 1.4; 
lam3 = 2.0;
lam4 = 1.3;

A = lam1*pts(1,:)';
B = lam2*pts(2,:)';
C = lam3*pts(4,:)';
D = lam4*pts(3,:)';
R = 1.025*sum(pts,1)'/8;


x = A(1) + ksi.*(B(1) - A(1)) + eta.*(D(1) - A(1)) + ksi.*eta.*((C(1)-A(1)) - (B(1)-A(1)) - (D(1) - A(1)));
y = A(2) + ksi.*(B(2) - A(2)) + eta.*(D(2) - A(2)) + ksi.*eta.*((C(2)-A(2)) - (B(2)-A(2)) - (D(2) - A(2)));
z = A(3) + ksi.*(B(3) - A(3)) + eta.*(D(3) - A(3)) + ksi.*eta.*((C(3)-A(3)) - (B(3)-A(3)) - (D(3) - A(3)));

normals = zeros(4,3);
ptss = [R';A';B';C';D'];

k = [1,2,3;1,3,4;1,4,5;1,5,2]; %back left front right
for i = 1:4
    tri = ptss(k(i,:),:);
    n1 = cross(tri(2,:) - tri(1,:), tri(3,:) - tri(1,:));
    tri = ptss(k(mod(i,4)+1,:),:);
    n2 = cross(tri(2,:) - tri(1,:), tri(3,:) - tri(1,:));
    n = (n1 + n2)/2;
    normals(i,:) = n / norm(n);
   
end
ids = [2,3,4,1];
normals =  -normals(ids,:);
original = ratio_vol(A,B,C,D,R);


figure 
axis off
hold on
view(-177,37)
trisurf(k,ptss(:,1),ptss(:,2),ptss(:,3),'FaceColor','green','FaceAlpha',0.5);
surf(x,y,z,'FaceAlpha',0.5);
str = {'R','A','B','C','D'};
for i = 1:5
    pt_i = ptss(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','r');
    text(pt_i(1),pt_i(2),pt_i(3),str{i}, 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
end
%exportgraphics(gca,'foto1.jpeg','Resolution',1000)

new_pts = zeros(5,3);
ksi_midle = [0.5,1,0.5,0,0.5];
eta_midle = [0,0.5,1.0,0.5,0.5];
x = A(1) + ksi_midle.*(B(1) - A(1)) + eta_midle.*(D(1) - A(1)) + ksi_midle.*eta_midle.*((C(1)-A(1)) - (B(1)-A(1)) - (D(1) - A(1)));
y = A(2) + ksi_midle.*(B(2) - A(2)) + eta_midle.*(D(2) - A(2)) + ksi_midle.*eta_midle.*((C(2)-A(2)) - (B(2)-A(2)) - (D(2) - A(2)));
z = A(3) + ksi_midle.*(B(3) - A(3)) + eta_midle.*(D(3) - A(3)) + ksi_midle.*eta_midle.*((C(3)-A(3)) - (B(3)-A(3)) - (D(3) - A(3)));

line([x(1),x(3)],[y(1),y(3)],[z(1),z(3)],'LineStyle','--');
line([x(2),x(4)],[y(2),y(4)],[z(2),z(4)],'LineStyle','--');
for i = 1:5
    new_pts(i,:) = [x(i),y(i),z(i)];
    plot3(x(i),y(i),z(i),'o','MarkerFaceColor','b','MarkerSize',5);
end
exportgraphics(gca,'foto2.jpeg','Resolution',1000)

ids = [3,4,1,2,5];
new_pts = new_pts(ids,:);
val1 = ratio_vol(new_pts(5,:)',new_pts(4,:)',C,new_pts(1,:)',R);
val2 = ratio_vol(new_pts(2,:)',new_pts(5,:)',new_pts(1,:)', D,R);
val3 = ratio_vol(A,new_pts(3,:)',new_pts(5,:)',new_pts(2,:)',R);
val4 = ratio_vol(new_pts(3,:)',B,new_pts(4,:)',new_pts(5,:)',R);
media = (val1 + val2 + val3 + val4)/4;
res = original/media;
hold off
figure
axis off
hold on
view(-177,37)
lam = 3;


draw_pyramid([new_pts(5,:); new_pts(4,:); C'; new_pts(1,:);R'] + lam*normals(1,:) );
draw_pyramid([new_pts(2,:); new_pts(5,:);new_pts(1,:); D' ;R'] + lam*normals(2,:));
draw_pyramid([A' ; new_pts(3,:); new_pts(5,:);new_pts(2,:);R'] + lam*normals(3,:) );
draw_pyramid([new_pts(3,:); B'; new_pts(4,:); new_pts(5,:);R'] + lam*normals(4,:));
exportgraphics(gca,'foto3.jpeg','Resolution',1000)

hold off


end
function out = ratio_vol(A,B,C,D,R)
volume_convexo = (1/6)*dot(A-R,cross((B-D),(C-A)));
volume_circ = (1/12)*dot(C-A,cross((D-A),(B-A)));
out =  (volume_circ)/volume_convexo;
out = abs(out);
end
function draw_pyramid(pts)
ksi = linspace(0,1,10);
eta = linspace(0,1,10);
[ksi,eta] = meshgrid(ksi,eta);
A = pts(1,:);
B = pts(2,:);
C = pts(3,:);
D = pts(4,:);
R = pts(5,:);
k = [1,2,3;1,3,4;1,4,5;1,5,2];
ptss = [R;A;B;C;D];

x = A(1) + ksi.*(B(1) - A(1)) + eta.*(D(1) - A(1)) + ksi.*eta.*((C(1)-A(1)) - (B(1)-A(1)) - (D(1) - A(1)));
y = A(2) + ksi.*(B(2) - A(2)) + eta.*(D(2) - A(2)) + ksi.*eta.*((C(2)-A(2)) - (B(2)-A(2)) - (D(2) - A(2)));
z = A(3) + ksi.*(B(3) - A(3)) + eta.*(D(3) - A(3)) + ksi.*eta.*((C(3)-A(3)) - (B(3)-A(3)) - (D(3) - A(3)));

trisurf(k,ptss(:,1),ptss(:,2),ptss(:,3),'FaceColor','green','FaceAlpha',0.5);
surf(x,y,z,'FaceAlpha',0.5);
end
function projected_coords = convex_study(pts, planes)

% pts = [0, 0, 0;10, 0, 10;0, 10, 10;10, 10, 0;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
% planes = [0.00000000  0.00000000  , 1.00000000  , -5.0;...
%     0.00000000  0.00000000  , -1.00000000 , 0;...
%     1.00000000  0.00000000  , 0.00000000  , -10;...
%     -1.00000000 0.00000000  , 0.00000000  , 0;...
%     0.00000000  -1.00000000 , 0.00000000  , 0;...
%     0.00000000  1.00000000  , 0.00000000  , -10];
% pts = [-2946.50000 , 2005.00000 , -3199.43994 ; ...
% -3034.50000 , 2055.00000 , -3220.27002 ; ...
% -3082.50000 , 1965.00000 , -3216.87012 ; ...
% -2994.50000 , 1925.00000 , -3194.76001 ; ...
% -2946.50000 , 2005.00000 , -3206.43994 ; ...
% -3034.50000 , 2055.00000 , -3227.27002 ; ...
% -3082.50000 , 1965.00000 , -3223.87012 ; ...
% -2994.50000 , 1925.00000 , -3201.76001 ];
% 
% planes = [-0.166330293 , 0.141016334  ,  0.975934744   , 2348.96509  ; ...
% 0.166330293  , -0.141016334 ,  -0.975934744  , -2355.79663  ; ...
% -0.882352948 , 0.470588237  , -0.00000000    , -3644.55884  ; ... 
% 0.857492924  , -0.514495790 ,  0.00000000    , 3558.16699  ; ... 
% 0.494009405  , 0.869456530  , -0.00000000    , -287.661621  ; ...  
% -0.413802952 , -0.910366535 , 0.00000000     , 513.322754  ];
% 
% pts = [-567.500000 , 2185.00000 , -3180.59009 ; ...
% -655.500000 , 2235.00000 , -3193.79004 ; ...
% -615.500000 , 2105.00000 , -3177.80005 ; ...
% -703.500000 , 2145.00000 , -3187.13989 ; ...
% -567.500000 , 2185.00000 , -3180.59009 ; ...
% -655.500000 , 2235.00000 , -3193.79004 ; ...
% -615.500000 , 2105.00000 , -3183.31006 ; ...
% -703.500000 , 2145.00000 , -3194.13989 ];
% 
% planes = [-0.0767133608 , 0.0976413265  , 0.992260635  ,  2899.79248 ; ...
% 0.112258554   , -0.0449573807 , -0.992661476 ,  -2995.77710  ; ...
% -0.882352948  , 0.470588237   , -0.00000000  ,  -1630.14709  ; ...
% 0.857492983   , -0.514495790  , 0.00000000   ,  1610.80054   ; ... 
% 0.413802952   , 0.910366535   , -0.00000000  ,  -1754.31775  ; ...
% -0.413802952  , -0.910366535  , 0.00000000   ,  1661.62585  ];
% 
% 
% pts = [441.140625 ,-1783.25000 , -3087.29004 ; ...
% 389.921875 ,-1695.75000 , -3084.12012 ; ...
% 354.578125 ,-1837.25000 , -3081.08008 ; ...
% 303.359375 ,-1749.75000 , -3072.65991 ; ...
% 441.140625 ,-1783.25000 , -3094.29004 ; ...
% 389.921875 ,-1695.75000 , -3085.14990 ; ...
% 354.578125 ,-1837.25000 , -3086.21997 ; ...
% 303.359375 ,-1749.75000 , -3072.65991 ];
% 
% planes = [0.104383603  , -0.00471933465 , 0.994525969    , 3017.22705 ;...
% -0.144522741 , 0.0435706489   , -0.988541782  , -2918.47021;...
% -0.529282808 , 0.848445415    , -0.00000000 , 1645.13025 ;...
% 0.529282928  , -0.848445415   , 0.00000000  , -1746.47852;...
% 0.863017738  , 0.505173564    , -0.00000000  , 520.138550 ;...
% -0.863017797 , -0.505173624   , 0.00000000  , -622.122986];



conec = [0,1,2;... %up
    1,3,2;...
    
    1,5,7;... %right
    7,3,1;...
    
    4,0,6;... %left
    0,2,6;...
    
    5,1,0;... %front
    0,4,5;...
    
    5,4,6;... %down
    5,6,7;...
    
    7,6,2;...%back
    7,2,3];

conec = [0, 3, 2; 0, 1, 3;
    4, 1, 0; 4, 5, 1;
    5, 3, 1; 5, 7, 3;...
    4, 2, 6; 4, 0, 2;
    7, 4, 6; 7, 5, 4;
    6, 3, 7; 6, 2, 3];
conec = conec + 1;

figure 
hold on 
axis (getBB(pts));
set(gca,'XColor','none','YColor','none','ZColor','none')
view(-85,57)
trisurf(conec,pts(:,1), pts(:,2), pts(:,3),'FaceAlpha',0.5,'FaceColor','r');
for i = 1:size(pts,1)
    pt_i = pts(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',4,'MarkerFaceColor','y');
    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i), 'VerticalAlignment','top','HorizontalAlignment','left');
end
%hold off

[k2,dualPoints,vol,bool,center] = findInterior2(planes,[0,0,0],false)
for i = 1:size(dualPoints,1)
    pt_i = dualPoints(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','cyan');
    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i), 'VerticalAlignment','top','HorizontalAlignment','left');
end
projected_coords = move(pts,planes);

view(-147,16)

exportgraphics(gca,strcat('proj_points_3','.png'),'Resolution',500);

return
j = 1;
for i = 3:size(planes,1)
    n = planes(i,1:3);
    d = planes(i,4);
    pt_i = pts(j,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',10,'MarkerFaceColor','y');
    drawPlan(n',d,'b');
    j = j + 1;
end
end
function new_coords = move(pts,planes)

% pts = [2999.82812 , -654.750000 , -3227.76001 ; ...
% 2911.82812 , -604.750000 , -3219.10010 ; ...
% 2863.82812 , -694.750000 , -3206.43994 ; ...
% 2951.82812 , -744.750000 , -3213.13989 ; ...
% 2999.82812 , -654.750000 , -3234.76001 ; ...
% 2911.82812 , -604.750000 , -3226.10010 ; ...
% 2863.82812 , -694.750000 , -3213.43994 ; ...
% 2951.82812 , -744.750000 , -3220.13989];
% 
% planes = [0.131468609 , 0.0796258524  ,0.988117337     ,2846.67236  ;...
% -0.131468624, -0.0796258673 ,-0.988117337    ,-2853.58887 ;...
% -0.882352948, 0.470588237   ,-0.00000000     ,2853.84839  ;...
% 0.882352948 , -0.470588237  ,0.00000000      ,-2955.02490 ;... 
% 0.494009405 , 0.869456530   ,-0.00000000     ,-912.666687 ;...  
% -0.494009405, -0.869456530  ,0.00000000      ,810.70306];

new_coords = zeros(8,3);
d_top = -planes(1,4);
d_bottom = -planes(2,4);
top = planes(1,1:3);
bottom = planes(2,1:3);
right = planes(3,1:3);
left = planes(4,1:3);
front = planes(5,1:3);
back = planes(6,1:3);
ids_pts = [1,2,4,3,5,6,8,7];
%pts = pts(ids_pts,:);
vert_planes = [left;front;right;back];
ids_planes = [4,5,3,6,1,2];
planes = planes(ids_planes,:);


for i = 1:4
p_i = pts(i,:);
p_j = pts(i+4,:);
edge = p_j - p_i;
if (dot(edge,edge) == 0)
    edge_dir = cross(vert_planes(i,:),vert_planes(mod(i,4)+1,:));
    p_i = p_j + edge_dir;
    p_j = p_j - edge_dir;
    edge = p_j - p_i;
end
n1 = vert_planes(i,:);
n2 = vert_planes(mod(i,4)+1,:) - dot(vert_planes(mod(i,4)+1,:),n1)*n1;
n2 = (1/dot(n2, vert_planes(mod(i,4)+1,:))) * n2;
n3_first = cross(n1,n2);
n3 = (1/dot(n3_first,top))*n3_first;

epslon = (5e-3);
d_i = dot(top,p_i) - d_top;
lam = -d_i / dot(top,edge);
new_coords(i,:) = p_i + lam*edge - ...
    epslon*(n1 + n2 + n3);

n3 = (1/dot(n3_first,bottom))*n3_first;

d_j = dot(bottom,p_j) - d_bottom;
lam = d_j / dot(bottom,edge);
new_coords(i + 4,:) = p_j - lam*edge - ...
    epslon*(n1 + n2 + n3);
end

for j = 1: 8
 center = new_coords(j,:);
for i = 1:size(planes,1)
    n = planes(i,1:3);
    d = planes(i,4);
    u = center + n * d;
    u = u / norm(u) ;
    angle = dot(u,n);
    if (angle >= 0.0)
         a = 1;
    end
    
end
end

% figure
% hold on
 axis (getBB(pts));
 set(gca,'XColor','none','YColor','none','ZColor','none')

% view(30,30);
n = planes(5,1:3);
d = planes(5,4)-0.2;
drawPlan(n',d,'b');
n = planes(6,1:3);
d = planes(6,4)-0.2;
drawPlan(n',d,'b');
view(-108,44)
%exportgraphics(gca,strcat('proj_points_1','.png'),'Resolution',500);

for i = 1:4
    pt_up_i = pts(i,:);
    pt_down_i = pts(i+4,:);
    pt_up_j = new_coords(i,:);
    pt_down_j = new_coords(i + 4,:);
    u_up = pt_up_j - pt_up_i;
    u_down = pt_down_j - pt_down_i;
    %plot3(pt_up_i(1),pt_up_i(2),pt_up_i(3),'o','MarkerSize',7,'MarkerFaceColor','b');
    %plot3(pt_down_i(1),pt_down_i(2),pt_down_i(3),'o','MarkerSize',7,'MarkerFaceColor','b');
    n = planes(i,1:3);
    d = planes(i,4)-0.2;
    %h1 = drawPlan(n',d,'b');
    n = planes(mod(i,4)+1,1:3);
    d = planes(mod(i,4)+1,4)-0.2;
    %h2 = drawPlan(n',d,'b');
    plot3(pt_up_j(1),pt_up_j(2),pt_up_j(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    plot3(pt_down_j(1),pt_down_j(2),pt_down_j(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    quiver3(pt_up_i(1),pt_up_i(2),pt_up_i(3),u_up(1),u_up(2),u_up(3),'g');
    quiver3(pt_down_i(1),pt_down_i(2),pt_down_i(3),u_down(1),u_down(2),u_down(3),'g');
    %delete([h1,h2]);
    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i), 'VerticalAlignment','top','HorizontalAlignment','left');
end
%exportgraphics(gca,strcat('proj_points_2','.png'),'Resolution',500);

end

function bb = getBB(coords)
margin = 1;
x_min = min(coords(:,1)) - margin;
x_max = max(coords(:,1)) + margin;
y_min = min(coords(:,2)) - margin;
y_max = max(coords(:,2)) + margin;
z_min = min(coords(:,3)) - margin;
z_max = max(coords(:,3)) + margin;
bb = [x_min,x_max,y_min,y_max,z_min,z_max];
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

function [k2,dualPoints,vol,bool,center] = findInterior2(M,trans,flag)
if nargin == 1
    trans = [0,0,0];
end
bool = false;
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
%plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','green');
if flag
    x = GetChebyshevCenter(M);
end
for i = 1:size(M,1)
    n = M(i,1:3);
    d = M(i,4);
    u = center + n * d;
    u = u / norm(u) ;
    angle = dot(u,n)
    if (angle >= 0.0)
        k2 = 0;
        dualPoints = [];
        vol = 0;
        return;
    end
    
end
bool = true;
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
% figure
% hold on
% view(30,30)
% trisurf(k1,coord_dual(:,1),coord_dual(:,2),coord_dual(:,3),'FaceColor','cyan');
% hold off
%volume = showDual(coord_dual,k1,center);
coord_points = coord_dual(unique(k1(:),'stable'),:);
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
    dualPoints(count,:) = (1/d)* normal + center + trans ;
    count = count + 1;
end
[k2,vol] = convhull(dualPoints);


% figure
% view(30,30);
% hold on
%trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
%hold off;
end