clc;
clear;
close all;

%BuildSubdivide();

%keyboard;
%find interior
%cell 76020
% showReal();
% keyboard;
% realVolume();
% keyboard;
% 
% 
% buildTetras()
% tetraPlan()
% keyboard
C_76020 = [3386.50000 , -1005.00000 , -3234.32007 ;...
    3299.50000 , -955.000000 , -3193.25000 ;...
    3338.50000 , -1085.00000 , -3217.97998 ;...
    3251.50000 , -1045.00000 , -3195.29004 ;...
    3386.50000 , -1005.00000 , -3234.32007 ;...
    3299.50000 , -955.000000 , -3193.25000 ;...
    3338.50000 , -1085.00000 , -3222.19995 ;...
    3251.50000 , -1045.00000 , -3195.29004 ];

C_36395 = [-459.500000 , -945.000000 , -3044.41992 ;...
    -547.500000 , -905.000000 , -3042.62988 ;...
    -507.500000 , -1035.00000 , -3048.44995 ;...
    -595.500000 , -985.000000 , -3051.25000 ;...
    -459.500000 , -945.000000 , -3044.41992 ;...
    -547.500000 , -905.000000 , -3049.62988 ;...
    -507.500000 , -1035.00000 , -3048.44995 ;...
    -595.500000 , -985.000000 , -3051.34009 ];

C_58569 = [824.500000  -1305.00000  -3070.21997 ;...
    736.500000  -1265.00000  -3077.97998 ;...
    776.500000  -1395.00000  -3058.84009 ;...
    688.500000  -1345.00000  -3076.92993 ;...
    824.500000  -1305.00000  -3070.21997 ;...
    736.500000  -1265.00000  -3084.97998 ;...
    776.500000  -1395.00000  -3058.84009 ;...
    688.500000  -1345.00000  -3076.92993 ];

C_32828 = [-380.500000 -1215.00000  -3038.76001 ;...
    -468.500000 -1175.00000  -3046.78003 ;...
    -428.500000 -1305.00000  -3040.09009 ;...
    -516.500000 -1255.00000  -3051.19995 ;...
    -380.500000 -1215.00000  -3038.76001 ;...
    -468.500000 -1175.00000  -3047.63989 ;...
    -428.500000 -1305.00000  -3040.09009 ;...
    -516.500000 -1255.00000  -3051.19995 ];

C_42085 = [-2698.50000 2675.00000  -3249.35010  ;...
    -2785.50000 2715.00000  -3246.71997  ;...
    -2746.50000 2585.00000  -3227.86011  ;...
    -2833.50000 2635.00000  -3235.37988  ;...
    -2698.50000 2675.00000  -3249.35010  ;...
    -2785.50000 2715.00000  -3253.71997  ;...
    -2746.50000 2585.00000  -3227.86011  ;...
    -2833.50000 2635.00000  -3235.37988  ];

C_64626 = [1465.50000  275.000000  -3119.10010  ;...
    1377.50000  325.000000  -3110.81006  ;...
    1417.50000  195.000000  -3122.96997  ;...
    1329.50000  235.000000  -3117.13989  ;...
    1465.50000  275.000000  -3119.10010  ;...
    1377.50000  325.000000  -3110.81006  ;...
    1417.50000  195.000000  -3125.08008  ;...
    1329.50000  235.000000  -3117.13989  ];

P_36395 = [-0.0338650532 -0.0564571023 , 0.997830451 ,2967.59155 ;...
    0.0490508452  0.00541325891 , -0.998781621,-3013.63428;...
    -0.857492924  0.514495790 	 , -0.00000000,-3.85870361;...
    0.857492924   -0.514495790  , 0.00000000 ,-92.1805115;...
    0.413802952   0.910366535 	 , -0.00000000,1050.43884 ;...
    -0.494009405  -0.869456530  , 0.00000000  ,-1150.59729]


P_58569 = [-0.0838609934 , 0.117119327  , 0.989570796 , 3257.88135 ; ...
    0.0975882709  , -0.166809827 , -0.981147885, -3310.08057; ...
    -0.857492924  , 0.514495790  , 0.00000000 , 1282.38062 ; ...
    0.857492924   , -0.514495790 , -0.00000000 , -1378.41992; ...
    0.413802952   , 0.910366535  , -0.00000000 , 846.847839 ; ...
    -0.413802952  , -0.910366535 , 0.00000000 , -939.539673];



P_32828 = [-0.0972447395 , 0.0204936918  , 0.995049536, 3010.88379  ; ...
    0.0990015790  , -0.0267710295 , -0.99472707 , -3017.09644; ...
    -0.857492924  , 0.514495730   , -0.00000000 , 202.797028 ; ...
    0.857492924   , -0.514495730  , 0.00000000  , -298.836243; ...
    0.413802952   , 0.910366476   , -0.00000000 , 1263.54724 ; ...
    -0.413802952  , -0.910366476  , 0.00000000  , -1356.23914]




P_42085 = [0.0547352023  , 0.156912580  , 0.986094534 , 2929.96313 ;...
    -0.0399984717 , -0.206201375 , -0.97769177 , -2732.94995;...
    -0.857492924  , 0.514495790  , -0.00000000 , -3785.40283;...
    0.857492924   , -0.514495790 , 0.00000000  , 3690.22070 ;...
    0.417733222   , 0.908569753  , -0.00000000 , -1303.17102;...
    -0.417733222  , -0.908569753 , 0.00000000  , 1210.43420 ]

P_64626 = [0.0382823236  , -0.0809300989 , 0.995984375 , 3072.32593 ;...
    -0.0425870381 , 0.0962222144  , -0.99444848 , -3065.98462;...
    -0.857492924  , 0.514495790   , -0.00000000 , 1019.13037 ;...
    0.857492924   , -0.514495790  , 0.00000000  , -1115.16956;...
    0.413802981   , 0.910366535   , -0.00000000 , -856.779053;...
    -0.413802981  , -0.910366535  , 0.00000000  , 764.087219 ]

coords = [2374.60840 -1284.10376 -2625.70093 ;...
    2309.30029 -1338.26208 -2646.57446 ;...
    2337.62988 -1343.49243 -2643.70752 ;...
    2240.23633 -1283.87805 -2644.42041 ;...
    2377.67554 -1286.07129 -2628.37793 ;...
    2312.39795 -1340.22888 -2649.28931 ;...
    2340.69189 -1345.44324 -2646.33862 ;...
    2241.75000 -1282.75000 -2646.97144 ];

coords = [-9775.53125 -39.7500000 -2965.29443 ;...
    -9855.40625 20.7500000  -2964.30103 ;...
    -9835.71875 -119.750000 -2965.60449 ;...
    -9915.59375 -59.2500000 -2963.61963 ;...
    -9775.53125 -39.7500000 -2972.03931 ;...
    -9855.40625 20.7500000  -2970.74390 ;...
    -9835.71875 -119.750000 -2971.65942 ;...
    -9915.59375 -59.2500000 -2969.81519 ];

coords = [-6541.54590 -3566.72437 -2745.17798 ;...
    -6638.37012 -3530.89526 -2747.98267 ;...
    -6543.19824 -3569.26758 -2745.17798 ;...
    -6640.13428 -3533.68799 -2747.98267 ;...
    -6542.00000 -3565.75000 -2746.87378 ;...
    -6638.68750 -3529.25000 -2749.93384 ;...
    -6541.99023 -3565.67358 -2746.87378 ;...
    -6638.68750 -3529.25000 -2749.93384 ];

test = [0.126878232   0.965406001   0.227800593     1574.66882  ; ...
    -0.129345447  -0.965394139  -0.226459816   -1566.80029 ; ...
    0.585910559   0.716896951   0.377845824   606.653870  ; ...
    0.380905688   -0.477089435  0.792020559   562.481689  ; ...
    -0.265553147  0.608599067   -0.747722387   -551.201050 ; ...
    -0.504498601  -0.813224494  -0.290081263   -680.516846 ; ...
    0.126878232   0.965406001   0.227800593     1574.66882  ; ...
    -0.129345447  -0.965394139  -0.226459816   -1566.80029 ; ...
    0.585910559   0.716896951   0.377845824   606.653870  ; ...
    0.380905688   -0.477089435  0.792020559   562.481689  ; ...
    -0.265553147  0.608599067   -0.747722387   -551.201050 ; ...
    -0.504498601  -0.813224494  -0.290081263   -680.516846 ];

test = [0.0129938750  -0.00745531591 0.999887884    3091.93555  ;...
    -0.0164629798 0.00420835195  -0.999855697  -3132.51465 ;...
    -0.799100101  0.601198018    -0.00000000   -7887.93066 ;...
    0.799100161   -0.601197958   0.00000000   7787.73096  ;...
    0.603785634   0.797146678    -0.00000000   5934.01172  ;...
    -0.603785694  -0.797146797   0.00000000   -6034.12451 ];

test = [-0.0234111417 0.0149990069  0.999613464   2644.46875  ;...
    -0.0301787108 0.00385624496 0.999537110   2561.92358  ;...
    -0.728837371  0.460405558   0.506776989   -1820.03381 ;...
    -0.285606414  -0.544970572  -0.788312137  -5977.04590 ;...
    0.288447887   0.812762558   0.506176591   6175.43896  ;...
    -0.113701768  -0.380523056  -0.917754948  -4621.58838 ];

P_1075506 = [0.333629042  0.399620771  0.853811920    2035.56738 ;...
    -0.333662361 -0.399652004 -0.853784382  -2035.14258;...
    1.00000000   0.00000000   0.00000000   3400.00000 ;...
    -1.00000000  -0.00000000  -0.00000000   -3500.00000;...
    -0.00000000  -1.00000000  -0.00000000   11600.0000 ;...
    0.00000000   1.00000000   0.00000000   -11700.0000];

C_1075506 = [-3500.00000 11600.0000 -6445.12012 ;...
    -3400.00000 11600.0000 -6485.47021 ;...
    -3500.00000 11700.0000 -6493.20020 ;...
    -3400.00000 11700.0000 -6531.00977 ;...
    -3500.00000 11600.0000 -6445.12012 ;...
    -3400.00000 11600.0000 -6485.47021 ;...
    -3500.00000 11700.0000 -6493.20020 ;...
    -3400.00000 11700.0000 -6531.02002 ];

P_61342 = [0.104383603  , -0.00471933465 , 0.994525969    , 3017.22705 ;...
-0.144522741 , 0.0435706489   , -0.988541782  , -2918.47021;...
-0.529282808 , 0.848445415    , -0.00000000 , 1645.13025 ;...
0.529282928  , -0.848445415   , 0.00000000  , -1746.47852;...
0.863017738  , 0.505173564    , -0.00000000  , 520.138550 ;...
-0.863017797 , -0.505173624   , 0.00000000  , -622.122986];



C_61342 = [441.140625 ,-1783.25000 , -3087.29004 ; ...
389.921875 ,-1695.75000 , -3084.12012 ; ...
354.578125 ,-1837.25000 , -3081.08008 ; ...
303.359375 ,-1749.75000 , -3072.65991 ; ...
441.140625 ,-1783.25000 , -3094.29004 ; ...
389.921875 ,-1695.75000 , -3085.14990 ; ...
354.578125 ,-1837.25000 , -3086.21997 ; ...
303.359375 ,-1749.75000 , -3072.65991 ];

%plotCell(C_1075506,P_1075506,1075506,'Proxy',30);

[k2,dualPoints,volume] = findInterior(P_61342,[0,0,0],false);

figure
view(30,30)
hold on
for i = 1:size(dualPoints,1)
    p = dualPoints(i,:);
    plot3(p(1),p(2),p(3),'o','MarkerFaceColor','yellow');
end
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
p = [303.35937500000000, -1749.7500000000000, -3073.7701157889633];
plot3(p(1),p(2),p(3),'o','MarkerFaceColor','blue');
Compare(C_61342,P_61342(1:6,:));
hold on
for i = 1:size(dualPoints,1)
    p = dualPoints(i,:);
    plot3(p(1),p(2),p(3),'o','MarkerFaceColor','yellow');
end
pi = [303.35937500000000, -1749.7500000000000, -3071.6603089014593];
pj = [303.35937500000000, -1749.7500000000000, -3073.6595153172907];
edge = pj - pi;
pj_k = pj + 1*edge;
pi_k = pi - 100*edge;
line([pi_k(1),pj_k(1)],[pi_k(2),pj_k(2)],[pi_k(3),pj_k(3)],'linewidth',3,'color','black');
plotCell(C_64626,P_64626,64626,'Pituba',5);

%planAndConvex(coords,test)

%plotCell(coords,test,1043882,'Prop',90);
%Compare(coords,test(1:6,:));
%Slides
case4 = [0, 0, 0;10, 0, 10;0, 10, 10;10, 10, 0;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
case4_trans = transform(case4);
cena = 1;
figure
set(gca,'XColor','none','YColor','none','ZColor','none')
hold on
axis ([0,20,0,20,0,15])
view(104,29)
plotSolid(case4,'red','yellow');
%exportgraphics(gca,strcat('Cena1','.png'),'Resolution',1000);
cena = cena+1;
plotSolid(case4_trans,'cyan','blue');
%exportgraphics(gca,strcat('Cena2','.png'),'Resolution',1000);
cena = cena + 1;
hold off
case4_planes = [0.00000000  0.00000000  , 1.00000000  , -5.0;...
    0.00000000  0.00000000  , -1.00000000 , 0;...
    1.00000000  0.00000000  , 0.00000000  , -10;...
    -1.00000000 0.00000000  , 0.00000000  , 0;...
    0.00000000  -1.00000000 , 0.00000000  , 0;...
    0.00000000  1.00000000  , 0.00000000  , -10];
%cena = ExportCompare(case4,case4_planes,cena,'r','y','g');
case4_trans_planes = [1.93572021e-08 ,  0.00000000  ,  1.00000000,-5.00000000;...
    0.00000000    ,  0.00000000  ,   -1.00000000,-0.00000000;...
    0.540301085   , 0.841471732  ,  0.00000000,-18.4146957;...
    -0.540301144  ,  -0.841471732,  0.00000000 ,8.41471672 ;...
    0.841471791   , -0.540301204 ,  0.00000000 ,-4.59696770;...
    -0.841471732  , 0.540301144  , 0.00000000 ,-5.40301132];
%cena = ExportCompare(case4_trans,case4_trans_planes,cena,'cyan','blue','magenta');

figure
set(gca,'XColor','none','YColor','none','ZColor','none')
hold on
axis ([0,20,0,20,0,15])
view(154,46)
for i = 1:size(case4_planes,1)
    n = case4_planes(i,1:3);
    d = case4_planes(i,4);
    drawPlan(n',d,'g')
end
for i = 1:size(case4_trans_planes,1)
    n = case4_trans_planes(i,1:3);
    d = case4_trans_planes(i,4);
    drawPlan(n',d,'magenta')
end
%exportgraphics(gca,strcat('Cena',num2str(cena),'.png'),'Resolution',1000);
cena = 18
cena = cena + 1;
hold off
[k2,dualPoints,volume,bool,center] = findInterior([case4_trans_planes;case4_planes],[0,0,0],false);
figure
set(gca,'XColor','none','YColor','none','ZColor','none')
hold on
axis ([0,20,0,20,0,15])
view(126,53)
for i = 1:size(case4_planes,1)
    n = case4_planes(i,1:3);
    d = case4_planes(i,4);
    drawPlan(n',d,'g')
end
for i = 1:size(case4_trans_planes,1)
    n = case4_trans_planes(i,1:3);
    d = case4_trans_planes(i,4);
    drawPlan(n',d,'magenta')
end
exportgraphics(gca,strcat('Cena',num2str(cena),'.png'),'Resolution',1000);
cena = cena + 1;
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
exportgraphics(gca,strcat('Cena',num2str(cena),'.png'),'Resolution',1000);
cena = cena + 1;
keyboard;




figure
planes =  [0.0700736716  , -0.0215827990,  0.997308314 ,3028.58496 ;...
    -0.0700736791 , 0.0215827953 ,  -0.997308314,-3035.56543;...
    -0.882352948  , 0.470588237  , -0.00000000 ,2172.79419 ;...
    0.882352948   ,-0.470588237  , 0.00000000 ,-2273.08838;...
    0.498283893   , 0.867013931  , -0.00000000,-2195.53857;...
    -0.498283893  , -0.867013931 ,  0.00000000 ,2093.58960 ];
A = [0.00000000  0.00000000  , 1.00000000  , -5.0;...
    0.00000000  0.00000000  , -1.00000000 , 0;...
    1.00000000  0.00000000  , 0.00000000  , -10;...
    -1.00000000 0.00000000  , 0.00000000  , 0;...
    0.00000000  -1.00000000 , 0.00000000  , 0;...
    0.00000000  1.00000000  , 0.00000000  , -10];
B = [1.93572021e-08 ,  0.00000000  ,  1.00000000,-5.00000000;...
    0.00000000    ,  0.00000000  ,   -1.00000000,-0.00000000;...
    0.540301085   , 0.841471732  ,  0.00000000,-18.4146957;...
    -0.540301144  ,  -0.841471732,  0.00000000 ,8.41471672 ;...
    0.841471791   , -0.540301204 ,  0.00000000 ,-4.59696770;...
    -0.841471732  , 0.540301144  , 0.00000000 ,-5.40301132];

hold on
[k2,dualPoints,volume] = findInterior([A;B],[0,0,0],false);
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');

keyboard

%% Halfplanes
conec = [0, 3, 2; 0, 1, 3;
    4, 1, 0; 4, 5, 1;
    5, 3, 1; 5, 7, 3;...
    4, 2, 6; 4, 0, 2;
    7, 4, 6; 7, 5, 4;
    6, 3, 7; 6, 2, 3];
conec_fitting = [0,1,3,2;...
    4,6,7,5;...
    5,7,3,1;... % um ponto
    4,0,2,6;...
    4,5,1,0;... % um ponto
    2,3,7,6];

case0 = [0, 0, 10;...
    10, 0, 10;...
    0, 10, 10;...
    10, 10, 10;...
    0, 0, 0;...
    10, 0, 0;...
    0, 10, 0;...
    10, 10, 0];
case1 = [0, 0, 10;0, 0, 10;0, 0, 10;0, 0, 10;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
case2 = [0, 0, 10;10, 0, 10;0, 10, 10;10, 10, 10;0, 0, 0;10, 0, 10;0, 10, 0;10, 10, 0];
case3 = [0, 0, 10;...
    10, 0, 10;...
    0, 10, 10;...
    10, 10, 10;...
    0, 0, 0;...
    10, 0, 10;...
    0, 10, 10;...
    10, 10, 0];
case4 = [0, 0, 0;...
    10, 0, 10;...
    0, 10, 10;...
    10, 10, 0;...
    0, 0, 0;...
    10, 0, 0;...
    0, 10, 0;...
    10, 10, 0];
case5 = [ 0, 0, 0;...  % 0
    0, 0, 0;... %% 1 [0,0,0] ou [10, 10, 0]
    0, 10, 10;... % 2
    10, 10, 0;... % 3
    0, 0, 0;...   % 4
    0, 0, 0;... %% 5 [0,0,0] ou [10, 10, 0]
    0, 10, 0;...  % 6
    10, 10, 0];   % 7

case_bruno =  [1276.25464 -2885.52686 -3043.70361;...
    1376.55835 -2885.95142 -3046.62451;...
    1279.13794 -2785.50562 -3049.44556;...
    1379.58325 -2786.01099 -3058.30469;...
    1276.25464 -2885.52686 -3043.70361;...
    1376.55835 -2885.95142 -3046.62451;...
    1279.13794 -2785.50562 -3049.44556;...
    1379.58325 -2786.01099 -3058.30469];

pts = {case1,case2,case3,case4,case5};
%% Volumes
real_volumes = [];
for i = 1:5
    figure
    hold on
    view(-172,33);
    volume = calculateVolume(pts{i},conec+1,i);
    real_volumes(end + 1) = volume;
    name = strcat('caseDegeneratedVolume',num2str(i));
    title(name);
    %exportgraphics(gca,strcat(name,'.png'),'Resolution',1000);
    hold off
end

for i = 1:5
    out = getVectorizedPoints(conec_fitting,pts{i});
end


A_trans = [-0.00000000, -0.00000000 , 1.00000000 , -10.0000000 ;...
    0.00000000 , 0.00000000  , -1.00000000           , -0.00000000 ;...
    0.382051021, 0.595010996 , 0.707106054           , -7.07106066 ;...
    -0.54030114, -0.841471732, 0.00000000            , -0.00000000 ;...
    0.841471732, -0.540301144, 0.00000000            , -0.00000000 ;...
    -0.59501105, 0.382050991 , 0.707106054           ,  -7.07106018 ];

B_trans = [0.00000000  , 0.00000000   , 1.00000000 ,  -10.0000000 ; ...
    0.367001027 , -0.367001027 , -0.854763448 			, -1.53310180 ; ...
    1.00000000  , 0.00000000   , 0.00000000   			, -10.0000000 ; ...
    -1.00000000 , 0.00000000   , 0.00000000   			, -0.00000000 ; ...
    0.00000000  , -1.00000000  , 0.00000000  			, -10.0000000 ; ...
    0.00000000  , 1.00000000   , 0.00000000  			, -0.00000000 ];

C_trans = [0.00000000  , 0.00000000  , 1.00000000 ,-10.0000000 ; ...
    0.00000000  , 0.00000000  , -1.00000000 ,5.00000000  ; ...
    1.00000000  , 0.00000000  , 0.00000000  ,-10.0000000 ; ...
    -1.00000000 , 0.00000000  , 0.00000000  ,-0.00000000 ; ...
    0.00000000  , -1.00000000 , 0.00000000  ,-10.0000000 ; ...
    0.00000000  , 1.00000000  , 0.00000000,-0.00000000 ];

D_trans = [1.93572021e-08 ,  0.00000000  ,  1.00000000,-5.00000000;...
    0.00000000    ,  0.00000000  ,   -1.00000000,-0.00000000;...
    0.540301085   , 0.841471732  ,  0.00000000,-18.4146957;...
    -0.540301144  ,  -0.841471732,  0.00000000 ,8.41471672 ;...
    0.841471791   , -0.540301204 ,  0.00000000 ,-4.59696770;...
    -0.841471732  , 0.540301144  , 0.00000000 ,-5.40301132];

E_trans = [0.797767580  , 0.173881054  , 0.577349424,-7.51230574;...
    0.00000000   , 0.00000000   , -1.00000000,-0.00000000;...
    -0.540301144  ,-0.841471732  , -0.00000000,8.41471672;...
    -0.540301144 , -0.841471732 , 0.00000000 ,8.41471672 ;...
    0.841471732  , -0.540301144 , -0.00000000,-4.59696770;...
    -0.841471732 , 0.540301144  , 0.00000,-5.40301132];


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

% E = [0.577350259 ,  -0.577350259 , 0.577350259 , -0.00000000 ;...
% 0.00000000  , 0.00000000    , -1.00000000 , -0.00000000 ;...
% 1.00000000  , -0.00000000   , -0.00000000 , -10.0000000 ;...
% -1.00000000 , 0.00000000   , 0.00000000  , -0.00000000 ;...
% -0.00000000 , -1.00000000   , -0.00000000 , 10.0000000 ;...
% 0.00000000  , 1.00000000    , 0.00000000  , -10.0000000 ];

normals = {A,B,C,D,E};
normals_trans = {[A;A_trans], [B;B_trans], [C;C_trans], [D;D_trans], [E;E_trans]};

% plotCellAndPlan(case5,conec + 1, E,conec_fitting);
%% Translation
pts_trans = {};
for i = 4:5
    pt_i = pts{i};
    pt_other = transform(pt_i);
    if (i == 2 || i == 3)
        pt_other = trans(pt_i);
    end
    %pts_trans(end+1) = pt_other;
    out = getVectorizedPoints(conec_fitting,pt_other);
    figure
    hold on
    view(149,39);
    plotSolid(pt_i,conec + ones(12,3),'cyan','blue');
    plotSolid(pt_other,conec + ones(12,3),'red','yellow');
    if i == 4
        normals = normals_trans{i};
        normals_1 = normals(1:6,:);
        normals_2 = normals(7:end,:);
        [k2,dualPoints_1,volume,valid] = findInterior(normals_1,[0,0,0],false);
        trisurf(k2,dualPoints_1(:,1),dualPoints_1(:,2),dualPoints_1(:,3),'FaceColor','cyan','FaceAlpha',0.2);
        [k2,dualPoints_2,volume,valid] = findInterior(normals_2,[0,0,0],false);
        trisurf(k2,dualPoints_2(:,1),dualPoints_2(:,2),dualPoints_2(:,3),'FaceColor','red','FaceAlpha',0.2);
    end
    [k2,dualPoints,volume,valid,center] = findInterior(normals_trans{i},[0,0,0],true);
    if valid
        trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','yellow','FaceAlpha',1);
    end
    name = strcat('TestsSolutions',num2str(i));
    title(name);
    %exportgraphics(gca,strcat(name,'.png'),'Resolution',1000);
    hold off
end


keyboard;

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
    trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
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
function plotSolid(pt,color,color_point)
conec = [0, 3, 2; 0, 1, 3;
    4, 1, 0; 4, 5, 1;
    5, 3, 1; 5, 7, 3;...
    4, 2, 6; 4, 0, 2;
    7, 4, 6; 7, 5, 4;
    6, 3, 7; 6, 2, 3];
conec = conec + 1;
for i = 1:size(pt,1)
    pt_i = pt(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',10,'MarkerFaceColor',color_point);
    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i), 'VerticalAlignment','top','HorizontalAlignment','left');
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
function plotPlane(A,pts)
for i = 1:size(A,1)
    n = A(i,1:3);
    d = A(i,4);
    pi = -n*d;
    value = realmax;
    n_pts = size(pts,1);
    p_end = zeros(1,3);
    pps = [];
    d = [];
    for j = 1:n_pts
        pj = pts(j,:) - pi;
        p_end_j = pj - dot(n,pj)*n + pi;
        p_end = p_end + p_end_j;
        pps(end+1, 1:3) = p_end_j;
        d(end + 1) = abs(dot(pj,n));
    end
    [~,ids] = sort(d);
    pps = pps(ids(1:4),:);
    m = mean(pps);
    lam = 1.4;
    xp1 = lam*(pps(1,:) - m) + m;
    value = 0;
    for j = 2:4
        d = norm(pps(j,:) - xp1);
        if d > value
            id = j;
            value = d;
        end
    end
    ids = 2:4;
    ids = ids(~ismember(ids,id));
    xp2 = lam*(pps(id,:) - m) + m;
    yp1 = lam*(pps(ids(1),:) - m ) + m;
    yp2 = lam*(pps(ids(2),:) - m) + m;
    %drawPlan(n',d,'r',100,p_end');
    h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
    [xp1(2),yp1(2),xp2(2),yp2(2)],...
    [xp1(3),yp1(3),xp2(3),yp2(3)],'r');
    set(h1, 'facealpha',0.3);
end
end
function plotPlaness(planes)
for i = 1:size(planes,1)
    n = planes(i,1:3);
    d = planes(i,4);
    drawPlan(n',d,'r');
end
end

function [k2,dualPoints,vol,bool,center] = findInterior(M,trans,flag)
if nargin == 1
    trans = [0,0,0];
    flag = false;
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
plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','green');
if flag
    x = GetChebyshevCenter(M);
end
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
figure
hold on
view(30,30)
trisurf(k1,coord_dual(:,1),coord_dual(:,2),coord_dual(:,3),'FaceColor','cyan');
hold off
volume = showDual(coord_dual,k1,center);
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
    %showTetra([pts;centroide]);
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
function plotCellAndPlan(cell,conec,normals,conec_fitting)
figure
hold on
view(86,44);
axis ([-1,11,-1,11,-1,11]);
[c,t] = getCellCentroide(cell,conec);
%[k2,dualPoints_1,volume,valid] = findInterior(normals,[0,0,0]);
%trisurf(k2,dualPoints_1(:,1),dualPoints_1(:,2),dualPoints_1(:,3),'FaceColor','red','FaceAlpha',0.2);
check = getVectorizedPoints(conec_fitting,cell);
[c,t] = getCellCentroide(cell,conec);
scale_cell = transpose(0.99*eye(3) * transpose(cell - c))+ c;
plotSolid(scale_cell,conec,'cyan','blue');
plotPlane(normals);
hold off
% figure
% hold on
% view(86,44);
% axis ([-1,11,-1,11,-1,11]);
% case1_scale = transpose(0.95*eye(3) * transpose(case1 + [0.1,0.1,0.1]));
% plotSolid(case1,conec + 1,'cyan','blue');
% plotPlane(A)
% %exportgraphics(gca,'halfplanes.png','Resolution',1000);
% hold off
end
function [centroide,triangles] = getCellCentroide(cell,conec)
n = size(conec,1);
triangles = [];
centroide = zeros(1,3);
count = 0;
for i = 1:n
    ids = conec(i,:);
    pts = cell(ids,:);
    u = pts(1,:) - pts(2,:);
    v = pts(1,:) - pts(3,:);
    area = 0.5*norm(cross(u,v));
    if area > 0
        triangles(end+1,1:3) = ids - 1;
        c = (pts(1,:) + pts(2,:) + pts(3,:))/3;
        centroide = centroide + c;
        count = count + 1;
    end
end
centroide = centroide / count;
end
function volume_other = showDual(cell,conec,center)

volume = 0.0;

volumes = [];
volume_other = 0;
volume_others = [];
ids = [6,0,5,1,8,2,11];

ids = ids + 1;

m = size(cell,1);

ids_i = 1:m;

bool = ~ismember(ids_i,ids);

ordem = [ids,ids_i(bool)];
%cell = [cell(ids,:);cell(bool,:)];

order = unique(conec(:),'stable');
figure
hold on
axis equal
view(30,30)

%plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','green');
areas = [];
for i = order'
    % find star
    star = [];
    for j = 1:size(conec,1)
        tri = conec(j,:);
        if ismember(i,tri)
            star(end+1,1:3) = tri;
        end
    end
    % put in cylick order
    star = setOrder(star,i);
    % find center of plan
    center = cell(i,:);
    d = norm(center);
    h = 1 / d;
    p0 = (center) * (h * h);
    plot3(p0(1),p0(2),p0(3),'o','MarkerSize',10);
    n = size(star,1)
    %teste
    volume_cell = 0;
    p_first = getDualPoint(cell(star(1,:),:));
    for l = 2:n-1
        pi = p_first;
        pj = getDualPoint(cell(star(l,:),:));
        pk = getDualPoint(cell(star(l+1,:),:));
        volume_cell = volume_cell + dot(pi, cross(pj,pk))
        pps = [pi;pj;pk];
        trisurf([1,2,3],pps(:,1),pps(:,2),pps(:,3));
    end
    volume_others(end+1) = volume_cell;
    volume_other = volume_other + volume_cell;
    area = 0;
    for k = 1:n
        pi = getDualPoint(cell(star(k,:),:));
        pj = getDualPoint(cell(star(mod(k,n)+1,:),:));
        u = pi - p0;
        v = pj - p0;
        normal = cross(u,v);
        area = area + norm(normal)*0.5;
        pt_triangles = [p0;pi;pj];
        trisurf([1,2,3],pt_triangles(:,1),pt_triangles(:,2),pt_triangles(:,3));
        
    end
    areas(end+1) = area;
    volumes(end+1) = area*(h/3);
    volume = volume + area*(h/3);
    
end
hold off
volume_other = volume_other/6;
end
function out = getDualPoint(tri)
p1 = tri(1,:);
p2 = tri(2,:);
p3 = tri(3,:);
normal = cross(p2 - p1, p3 - p1);
normal = normal / norm(normal);
d = dot(normal,p1);
out = (1/d)* normal ;
end
function star = setOrder(tris,i)
star = [];
first_tri = tris(1,:);
star(end+1,1:3) = first_tri;
[~,pos] = ismember(i,first_tri);
next = mod(pos,3)+1;
bef = mod(next,3)+1;
ids_2 = first_tri(next);
ids_3 = first_tri(bef);
fixed = ids_2;
next = ids_3;
while next ~= fixed
    [row,col] = find(tris == i);
    next_array = mod(col,3) + 1;
    [next_tri,~] = find(diag(tris(row,next_array)) == next);
    selected_row = row(next_tri);
    tri_i = tris(selected_row,:);
    next = tri_i(mod(next_array(next_tri) , 3) + 1);
    star(end+1,1:3) = tri_i;
end

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
function plotCell(coords,planes,id,filename,erro)
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

conec = conec + 1;
volume = 0 ;
figure
view(30,30)
hold on
if ~isempty(coords)
    for i = 1:12
        conec_i = conec(i,:);
        pi = coords(conec_i(1),:);
        pj = coords(conec_i(2),:);
        pk = coords(conec_i(3),:);
        volume_i  = dot(pi, cross(pj,pk))/6
        pts = [pi;pj;pk];
        volume = volume + abs(volume_i);
        trisurf([1,2,3],pts(:,1),pts(:,2),pts(:,3),'FaceColor','cyan');
        
    end
end
hold off
figure
hold on
[k2,dualPoints,volume] = findInterior(planes,[0,0,0],false);
volume
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
hold off

fig = figure
hold on
view(70,13);
axis (getBB(coords))
new_filename = strcat(filename,',','erro:',num2str(erro),'%','id:');
name = strcat(new_filename,num2str(id));

title(name,'color','black');
set(gca,'XColor','none','YColor','none','ZColor','none')
%set(gca,'title',title_text);
set(gcf,'color','white');
trisurf(conec,coords(:,1),coords(:,2),coords(:,3),'FaceColor','cyan','FaceAlpha',0.5);
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red','FaceAlpha',0.5);
legend('exato','convexo');
%text(10,10,name,'color','black');
exportgraphics(gca,strcat(filename,num2str(id),'.png'),'Resolution',1000);
hold off
end
function Compare(coords,planes)
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
conec = conec + 1;
figure
view(70,-13)
hold on
set(gca,'XColor','none','YColor','none','ZColor','none')
axis (getBB(coords))
[c,t] = getCellCentroide(coords,conec);
scale_cell = transpose(0.95*eye(3) * transpose(coords - c ))+ c - [0,0,0.5];
%trisurf(conec,scale_cell(:,1),scale_cell(:,2),scale_cell(:,3),'FaceColor','cyan');
trisurf(conec,coords(:,1),coords(:,2),coords(:,3),'FaceColor','cyan');

for i = 1:2
    n = planes(i,1:3);
    d = planes(i,4);
    drawPlan(n',d,'r')
end
%exportgraphics(gca,strcat('Pituba_planes','.png'),'Resolution',1000);

A = planes(:,1:3);
A(:,4) = 1.0;
b = -planes(:,4);
n = size(A,2);
m = size(A,1);
c = zeros(n,1);
c(end) = -1;
tic
center = linprog(c,A,b);
r = center(4,1);
center = linprog(c,A,b);
r = center(4,1);
center = center(1:3)';
plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','green');
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
function planAndConvex(coords,planes)
figure
view(30,30)
[k2,dualPoints,volume,bool,center] = findInterior(planes,[0,0,0],false);

bb = getBB(dualPoints);
margem = [-100 100 -100 100 -100 100];
bb = bb + margem;
hold on
axis (bb)
plot3(center(1),center(2),center(3),'o','MarkerSize',5,'MarkerFaceColor','blue');
trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','cyan');

for i = 1:size(planes)
    n = planes(i,1:3);
    d = planes(i,4);
    drawPlan(n',d,'r')
end
hold off
end
function cena = ExportCompare(coords,planes,cena,color1,color2,color3)
figure
hold on
axis ([0,20,0,20,0,15])
view(163,42)
set(gca,'XColor','none','YColor','none','ZColor','none')
plotSolid(coords,color1,color2);
exportgraphics(gca,strcat('Cena',num2str(cena),'.png'),'Resolution',1000);
cena = cena + 1;
for i = 1:size(planes,1)
    n = planes(i,1:3);
    d = planes(i,4);
    drawPlan(n',d,color3)
    exportgraphics(gca,strcat('Cena',num2str(cena),'.png'),'Resolution',1000);
    cena = cena + 1;
end
hold off
end
function tetraPlan()
C_146 = [3005.50000 , 805.000000 , -3230.40991 ;...
    2918.50000 , 855.000000 , -3223.41992 ;...
    2957.50000 , 715.000000 , -3229.18994 ;...
    2870.50000 , 765.000000 , -3221.79004 ;...
    3005.50000 , 805.000000 , -3237.40991 ;...
    2918.50000 , 855.000000 , -3230.41992 ;...
    2957.50000 , 715.000000 , -3236.18994 ;...
    2870.50000 , 765.000000 , -3228.79004 ];
C_227 = [2782.50000, 815.000000 ,-3217.43994 ;...
    2694.50000, 855.000000 ,-3212.27002 ;...
    2734.50000, 725.000000 ,-3204.87012 ;...
    2646.50000, 775.000000 ,-3196.41992 ;...
    2782.50000, 815.000000 ,-3224.43994 ;...
    2694.50000, 855.000000 ,-3219.27002 ;...
    2734.50000, 725.000000 ,-3211.87012 ;...
    2646.50000, 775.000000 ,-3203.41992 ];
C_52436 = [3005.50000 , 805.000000 , -3321.67993;...
    2918.50000 , 855.000000 , -3303.70996;...
    2957.50000 , 715.000000 , -3327.18994;...
    2870.50000 , 765.000000 , -3304.60010;...
    3005.50000 , 805.000000 , -3321.67993;...
    2918.50000 , 855.000000 , -3303.70996;...
    2957.50000 , 715.000000 , -3327.64990;...
    2870.50000 , 765.000000 , -3304.60010];
case0 = [0, 0, 10;...
    10, 0, 10;...
    0, 10, 10;...
    10, 10, 10;...
    0, 0, 10;...
    10, 0, 10;...
    0, 10, 0;...
    10, 10, 10];
P_146 = [0.0700736716  , -0.0215827990 , 0.997308314  ,     3028.58496  ;...
    -0.0700736791 , 0.0215827953  , -0.997308314 ,  -3035.56543 ;...
    -0.882352948  , 0.470588237   , -0.00000000  ,  2172.79419  ;...
    0.882352948   , -0.470588237  , 0.00000000   , -2273.08838 ;...
    0.498283893   , 0.867013931   , -0.00000000  ,  -2195.53857 ;...
    -0.498283893  , -0.867013931  , 0.00000000   , 2093.58960  ];
P_227 = [-0.124739498 , -0.0957080945 , -0.987562597 ,  -2753.39990 ;...
    0.124739513  , 0.0957081020  , 0.987562656  ,  2760.31299  ;...
    0.857492924  , -0.514495790  , 0.00000000   , -1870.62073 ;...
    -0.882352948 , 0.470588237   , -0.00000000  , 2071.61768  ;...
    -0.413802952 , -0.910366535  , 0.00000000   , 1893.35547  ;...
    0.494009405  , 0.869456530   , -0.00000000  , -1981.22473 ];
%[k2,dualPoints,volume] = findInterior(P_227,[0,0,0],false);
P_52436 = [-0.159422010 , 0.119853385  , -0.979908168 , -2873.40796 ;...
    0.160205841  , -0.122752592 , 0.979421139  , 2871.87915  ;...
    0.882352948  , -0.470588237 , 0.00000000   , -2172.79419 ;...
    -0.882352948 , 0.470588237  , -0.00000000  , 2273.08838  ;...
    -0.498283893 , -0.867013991 , 0.00000000   , 2195.53857  ;...
    0.498283893  , 0.867013991  , -0.00000000  , -2093.58960 ];
[k2,dualPoints,volume,bool,center] = findInterior(-P_52436,[0,0,0],false);

bb = getBB(C_52436);
L = [bb(2)-bb(1),bb(4)-bb(3),bb(6)-bb(5)];
%Compare(C_52436,P_52436(1:6,:));
tetras = plotTetras(case0,1);
ids_tetra = [1, 3, 2;...
    0, 2, 3;...
    0, 3, 1;...
    0, 1, 2];
ids_tetra = ids_tetra + 1;
for i = 1:6
    dists = clipTetraPlan(tetras{i},P_52436);
    volume = calculateVolume(tetras{i},ids_tetra,-1);
end
end
function tetras = plotTetras(coords,n)
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
conec = conec + 1;
figure
view(70,-13)
hold on
set(gca,'XColor','none','YColor','none','ZColor','none')
axis (getBB(coords))
[c,t] = getCellCentroide(coords,conec);
scale_cell = transpose(0.95*eye(3) * transpose(coords - c ))+ c - [0,0,0.5];
scale_cell = coords;
trisurf(conec,scale_cell(:,1),scale_cell(:,2),scale_cell(:,3),'FaceColor','cyan','FaceAlpha',0.5);
ids = [0, 1, 5, 7;...
    0, 5, 4, 7;...
    0, 4, 6, 7;...
    0, 6, 2, 7;...
    0, 2, 3, 7;...
    0, 3, 1, 7];
ids = ids + 1;
ids_tetra = [1, 3, 2;...
    0, 2, 3;...
    0, 3, 1;...
    0, 1, 2];
ids_tetra = ids_tetra + 1;
tetras = cell(6,1);
for i = 1:n
    id = ids(i,:);
    vertices = coords(id,:);
    tetras{i} = vertices;
    plot3(vertices(:,1),vertices(:,2),vertices(:,3),'o','MarkerFaceColor','yellow','MarkerSize',10);
    trisurf(ids_tetra,vertices(:,1),vertices(:,2),vertices(:,3),'FaceColor','red');
end

end
function dists = clipTetraPlan(vertices,planes)
dists = zeros(4,6);
for i = 1:4
    v = vertices(i,:);
    for j = 1:6
        n = planes(j,1:3);
        d = planes(j,4);
        u = v + n*d;
        angle = dot(n,u);
        dists(i,j) = angle;
    end
end
end
function buildTetras()
case4 = [0, 0, 0;10, 0, 10;0, 10, 10;10, 10, 0;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
case1 = [0, 0, 10;0, 0, 10;0, 0, 10;0, 0, 10;0, 0, 0;10, 0, 0;0, 10, 0;10, 10, 0];
% case4 = [3005.50000 , 805.000000 , -3321.67993;...
%     2918.50000 , 855.000000 , -3303.70996;...
%     2957.50000 , 715.000000 , -3327.18994;...
%     2870.50000 , 765.000000 , -3304.60010;...
%     3005.50000 , 805.000000 , -3321.67993;...
%     2918.50000 , 855.000000 , -3303.70996;...
%     2957.50000 , 715.000000 , -3327.64990;...
%     2870.50000 , 765.000000 , -3304.60010];
case4 = [809.500000 , -505.000000 , -3107.82007;...
721.500000 , -455.000000 , -3109.54004 ;...
761.500000 , -585.000000 , -3096.38989 ;...
673.500000 , -545.000000 , -3102.27002 ;...
809.500000 , -505.000000 , -3108.11011 ;...
721.500000 , -455.000000 , -3109.54004 ;...
761.500000 , -585.000000 , -3096.38989 ;...
673.500000 , -545.000000 , -3102.27002 ];
case4 = [-2881.50000 , 2545.00000 , -3229.04004 ;...
-2969.50000 , 2595.00000 , -3240.97998 ;...
-2929.50000 , 2455.00000 , -3252.80005 ;...
-3017.50000 , 2505.00000 , -3235.60010 ;...
-2881.50000 , 2545.00000 , -3229.04004 ;...
-2969.50000 , 2595.00000 , -3246.05005 ;...
-2929.50000 , 2455.00000 , -3259.80005 ;...
-3017.50000 , 2505.00000 , -3242.60010 ];

%showCell(case4);
conec = [0, 3, 2; 0, 1, 3;
    4, 1, 0; 4, 5, 1;
    5, 3, 1; 5, 7, 3;...
    4, 2, 6; 4, 0, 2;
    7, 4, 6; 7, 5, 4;
    6, 3, 7; 6, 2, 3];
heights = [case4(1,:) - case4(5,:); case4(2,:) - case4(6,:);case4(3,:)-case4(7,:);case4(4,:) - case4(8,:)];
norms = diag(heights * heights');
out = 1;
conec = conec + 1;
colision_1 = 5;
colision_2 = 1;
tetras = [];
for i = 1:size(conec,1)
    ids = conec(i,:);
    if ismember(colision_1,ids) || ismember(colision_2,ids)
        continue
    end
    ids_l = ids - 1;
    %is degenrate?
    pts = case4(ids,:);
    u = pts(1,:) - pts(2,:);
    v = pts(1,:) - pts(3,:);
    n = cross(u,v);
    area = dot(n,n);
    %tetras(end+1,1:4) = [ids,5];
    if area > 0
        tetras(end+1,1:4) = [ids,5];
    end
end
%tetras = [];
%tetras(1:2,1:4) = [  8, 7, 6,out; 7,5,6,out];
%tetras(1:2,1:4) = [  8, 7, 5,out; 5,6,8,out];

ids_tetra = [1, 3, 2;...
    0, 2, 3;...
    0, 3, 1;...
    0, 1, 2];
ids_tetra = ids_tetra + 1;
hold on
view(30,30)
volume = 0;
for i = 1:size(tetras,1)
    vertices = case4(tetras(i,:),:);
    u = vertices(1,:) - vertices(2,:);
    v = vertices(1,:) - vertices(3,:);
    n = cross(u,v);
    norm_n = norm(n);
    n = n / norm_n;
    area = norm_n;
    d = vertices(4,:) -  vertices(1,:);
    centroide = vertices(4,:);
    h = abs(dot(n,d));
    volume = volume + h*area*(1/6);
    %showTetra([vertices(1:3,:);vertices(4,:)]);
    plot3(vertices(:,1),vertices(:,2),vertices(:,3),'o','MarkerFaceColor','yellow','MarkerSize',10);
    plot3(vertices(4,1),vertices(4,2),vertices(4,3),'o','MarkerFaceColor','black','MarkerSize',10);
    trisurf(ids_tetra,vertices(:,1),vertices(:,2),vertices(:,3));
end


end
function showCell(coords)
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

conec = conec + 1;
trisurf(conec,coords(:,1), coords(:,2), coords(:,3),'FaceColor','cyan');
end
function showHalfSpace(pts)
p1 = [pts(1,:);pts(2,:);pts(4,:); pts(3,:)];

p2 = [pts(5,:);pts(7,:);pts(8,:); pts(6,:)];

p3 = [pts(1,:);pts(3,:);pts(7,:); pts(5,:)];

p4 = [pts(2,:);pts(6,:);pts(8,:); pts(4,:)];

p5 = [pts(1,:);pts(5,:);pts(6,:); pts(2,:)];

p6 = [pts(3,:);pts(4,:);pts(8,:); pts(7,:)];

hold on
showCell(pts);
%NormalsFromTopsPlane(p1);
%NormalsFromTopsPlane(p2);
%NormalsFromVerticalPlanes(p3);
%NormalsFromVerticalPlanes(p4);
%NormalsFromVerticalPlanes(p5);
%NormalsFromVerticalPlanes(p6,true);
end
function NormalsFromVerticalPlanes(plane, flag)
if nargin == 1
    flag = false;
end
t1 = [1,2,4];
p_t1 = plane(t1,:);

n1 = cross(p_t1(1,:) - p_t1(2,:), p_t1(1,:) - p_t1(3,:));
n1 = 30*n1/norm(n1);
c_t1 = sum(plane)/4;
if flag 
    n1 = n1;
end
% if (dot(n1,p_t1(1,:)) < 0)
%     n1 = - n1;
% end
plot3(c_t1(1),c_t1(2),c_t1(3),'o','MarkerSize',4,'MarkerFaceColor','red');
quiver3(c_t1(1),c_t1(2),c_t1(3),n1(1),n1(2),n1(3),'color','black');

end
function [n1,n2] = NormalsFromTopsPlane(plane)
t1 = [1,2,4];
t2 = [3,4,2];
p_t1 = plane(t1,:);
p_t2 = plane(t2,:);

n1 = cross(p_t1(1,:) - p_t1(2,:), p_t1(1,:) - p_t1(3,:));
n2 = cross(p_t2(1,:) - p_t2(2,:), p_t2(1,:) - p_t2(3,:));
n1 = 10*n1/norm(n1);
n2 = 10*n2/norm(n2);

c_t1 = sum(p_t1)/3;
c_t2 = sum(p_t2)/3;

quiver3(c_t1(1),c_t1(2),c_t1(3),n1(1),n1(2),n1(3));
quiver3(c_t2(1),c_t2(2),c_t2(3),n2(1),n2(2),n2(3));

end
function showReal()
 pts = [809.500000 , -505.000000 , -3107.82007;...
721.500000 , -455.000000 , -3109.54004 ;...
761.500000 , -585.000000 , -3096.38989 ;...
673.500000 , -545.000000 , -3102.27002 ;...
809.500000 , -505.000000 , -3108.11011 ;...
721.500000 , -455.000000 , -3109.54004 ;...
761.500000 , -585.000000 , -3096.38989 ;...
673.500000 , -545.000000 , -3102.27002 ];
pts = [0, 0, 10;...
    10, 0, 10;...
    0, 10, 10;...
    10, 10, 10;...
    0, 0, 0;...
    10, 0, 0;...
    0, 10, 0;...
    10, 10, 0];
pts = pts + [0,0,1];
pts(1,:) = [pts(1,1),pts(1,2),1.65*pts(1,3)];
pts(2,:) = [pts(2,1),pts(2,2),1.25*pts(2,3)];
pts(3,:) = [pts(3,1),pts(3,2),1.15*pts(3,3)];
pts(4,:) = [pts(4,1),pts(4,2),0.95*pts(4,3)];

pts(5,:) = [pts(5,1),pts(5,2),0.75*pts(5,3)];
pts(6,:) = [pts(6,1),pts(6,2),1.95*pts(6,3)];
pts(7,:) = [pts(7,1),pts(7,2),1.65*pts(7,3)];
pts(8,:) = [pts(8,1),pts(8,2),0.8*pts(8,3)];


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
figure 
axis off
hold on
view(-18,-23)
%text
for i = 1:size(pts,1)
    pt_i = pts(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','r');
    text(pt_i(1),pt_i(2),pt_i(3),num2str(i-1), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
end
k = [1,2,3;1,3,4;1,4,5;1,5,2];
R = 1.025*sum(pts,1)/8;
%plot3(R(1),R(2),R(3),'o','MarkerSize',10,'MarkerFaceColor','b');
for i = 1:6
    coords = planes{i};
    p0 = coords{1}';
    p1 = coords{2}';
    p2 = coords{3}';
    p3 = coords{4}';
    
    x = p0(1) + ksi.*(p1(1) - p0(1)) + eta.*(p3(1) - p0(1)) + ksi.*eta.*( (p2(1)-p0(1)) - (p1(1)-p0(1)) - (p3(1) - p0(1)) );
    y = p0(2) + ksi.*(p1(2) - p0(2)) + eta.*(p3(2) - p0(2)) + ksi.*eta.*((p2(2)-p0(2)) - (p1(2)-p0(2)) - (p3(2) - p0(2)));
    z = p0(3) + ksi.*(p1(3) - p0(3)) + eta.*(p3(3) - p0(3)) + ksi.*eta.*((p2(3)-p0(3)) - (p1(3)-p0(3)) - (p3(3) - p0(3)));
    if i == 1 || i == 2
        %surf(x,y,z,'FaceAlpha',0.5,'EdgeAlpha',0.1, 'FaceColor','red');
    else
        surf(x,y,z,'FaceAlpha',0.5,'EdgeAlpha',0.1);
    end
    %k1 = [2,3,4;2,4,5];
    k1 = [2,3,5;3,4,5];
    ptss = [R;p0';p1';p2';p3'];
    if i == 1 || i == 2
        if i == 1
            k1 = [2,3,4;2,4,5];
            %k1 = [2,3,5;3,4,5];
        else
            %k1 = [2,3,4;2,4,5];
            k1 = [2,3,5;3,4,5];
        end
        trisurf(k1,ptss(:,1),ptss(:,2),ptss(:,3),'FaceColor','blue','FaceAlpha',0.5);
    end

    for i = 1:4
        pi = coords{i}';
        pj = coords{mod(i,4)+1}';
        line([pi(1),pj(1)],[pi(2),pj(2)],[pi(3),pj(3)],'color','black');
    end
    %pyramid
    pyramid = [R;p0';p1';p2';p3'];
    %trisurf(k,pyramid(:,1),pyramid(:,2),pyramid(:,3),'FaceColor','green','FaceAlpha',0.2);

end
exportgraphics(gca,'geral_planes_3.png','Resolution',1000)
hold off

c = 1.025*sum(pts,1)/8;
p0 = pts(1,:)';
p1 = pts(2,:)';
p2 = pts(4,:)';
p3 = 1.01*pts(3,:)';


x = p0(1) + ksi.*(p1(1) - p0(1)) + eta.*(p3(1) - p0(1)) + ksi.*eta.*((p2(1)-p0(1)) - (p1(1)-p0(1)) - (p3(1) - p0(1)));
y = p0(2) + ksi.*(p1(2) - p0(2)) + eta.*(p3(2) - p0(2)) + ksi.*eta.*((p2(2)-p0(2)) - (p1(2)-p0(2)) - (p3(2) - p0(2)));
z = p0(3) + ksi.*(p1(3) - p0(3)) + eta.*(p3(3) - p0(3)) + ksi.*eta.*((p2(3)-p0(3)) - (p1(3)-p0(3)) - (p3(3) - p0(3)));

k = [1,2,3;1,3,4;1,4,5;1,5,2];
k1 = [2,3,4;2,4,5];
ptss = [c;p0';p1';p2';p3'];
figure 
axis off
hold on
view(170,11)
trisurf(k,ptss(:,1),ptss(:,2),ptss(:,3),'FaceColor','green','FaceAlpha',0.5);
%trisurf(k1,pts(:,1),pts(:,2),pts(:,3),'FaceColor','blue','FaceAlpha',0.5);
surf(x,y,z,'FaceAlpha',0.5);
str = {'R','A','B','C','D'};
for i = 1:5
    pt_i = ptss(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','r');
    text(pt_i(1),pt_i(2),pt_i(3),str{i}, 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
end
%exportgraphics(gca,'pyramid.png','Resolution',1000)
hold off
figure
axis off
hold on
view(170,11)
trisurf(k,ptss(:,1),ptss(:,2),ptss(:,3),'FaceColor','green','FaceAlpha',0.5);
trisurf(k1,ptss(:,1),ptss(:,2),ptss(:,3),'FaceColor','blue','FaceAlpha',0.5);
for i = 1:5
    pt_i = ptss(i,:);
    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','r');
    text(pt_i(1),pt_i(2),pt_i(3),str{i}, 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
end
%exportgraphics(gca,'linear_pyramid.png','Resolution',1000)
hold off
figure
axis off
hold on
view(37,25)
for i = 1:size(pts,1)
    pt_i = pts(i,:);
    %plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','r');
    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i-1), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
end

%Linear aprox
normals = zeros(6,4);
for i = 1:6
    coords = planes{i};
    A = coords{1};
    B = coords{2};
    C = coords{3};
    D = coords{4};
    ptss = [A;B;C;D];
    n1 = cross((B-A),(C-A));
    n2 = cross((C-A),(D-A));
    n1 = n1/norm(n1);
    n2 = n2/norm(n2);
    n = (n1 + n2)/2;
    n =  n/norm(n);
    d = dot(n,A);
    normals(i,1:3) = n/norm(n);
    normals(i,4) = d;
    %trisurf([1,2,3;3,4,1],ptss(:,1),ptss(:,2),ptss(:,3),'FaceColor','blue','FaceAlpha',0.2);
end
normals(1,1:3) = -normals(1,1:3);
normals(4,1:3) = -normals(4,1:3);
normals(6,1:3) = -normals(6,1:3);


%exportgraphics(gca,'linear_cell.png','Resolution',1000)
planess = cell(3,1);
planess{1} = {pts(5,:),pts(7,:),pts(8,:), pts(6,:)};
planess{2} = {pts(2,:),pts(6,:),pts(8,:), pts(4,:)};
planess{3} = {pts(3,:),pts(4,:),pts(8,:), pts(7,:)};
apex = pts(1,:);
k = [1,2,4;2,3,4;3,1,4];
plot3(apex(1),apex(2),apex(3),'o','MarkerSize',15,'MarkerFaceColor','r');
for i = 1:3
    coords = planess{i};
    A = coords{1};
    B = coords{2};
    C = coords{3};
    D = coords{4};
    tetra1 = [C;B;A;apex];
    tetra2 = [D;C;A;apex];
    trisurf(k,tetra1(:,1),tetra1(:,2),tetra1(:,3),'FaceColor','red','FaceAlpha',0.2);
    trisurf(k,tetra2(:,1),tetra2(:,2),tetra2(:,3),'FaceColor','red','FaceAlpha',0.2);
    trisurf([1,2,3],tetra1(:,1),tetra1(:,2),tetra1(:,3),'FaceColor','blue','FaceAlpha',0.6);
    trisurf([1,2,3],tetra2(:,1),tetra2(:,2),tetra2(:,3),'FaceColor','blue','FaceAlpha',0.6);
end
%exportgraphics(gca,'linear_cell_pyramid.png','Resolution',1000)

hold off
%HalfPlanes
lims = getBB(pts);
figure
hold on
axis ([lims(1)-1,lims(2)+1,lims(3)-1,lims(4)+1,lims(5)-1,lims(6)+1])
axis off
plotPlane(normals);
view(82,21)
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
    for i = 1:4
        pi = coords{i}';
        pj = coords{mod(i,4)+1}';
        line([pi(1),pj(1)],[pi(2),pj(2)],[pi(3),pj(3)],'color','black');
    end
end
exportgraphics(gca,'half.png','Resolution',1000)
hold off



end
function [volume_total,volume_concave,ratio, volume_convex] = realVolume(pts)
planes{1} = {pts(1,:),pts(2,:),pts(4,:), pts(3,:)};
planes{2} = {pts(5,:),pts(7,:),pts(8,:), pts(6,:)};
planes{3} = {pts(1,:),pts(3,:),pts(7,:), pts(5,:)};
planes{4} = {pts(2,:),pts(6,:),pts(8,:), pts(4,:)};
planes{5} = {pts(1,:),pts(5,:),pts(6,:), pts(2,:)};
planes{6} = {pts(3,:),pts(4,:),pts(8,:), pts(7,:)};
R = sum(pts,1)/8;
R = R';
volume_total = 0.0;
volume_concave = 0.0;
volume_convex = 0.0;
for i = 1:6
    coords = planes{i};
    A = coords{1}';
    B = coords{2}';
    C = coords{3}';
    D = coords{4}';
    volume = (1/6)*dot(A-R,cross((B-D),(C-A))) + (1/12)*dot(C-A,cross((D-A),(B-A)));
    volume_total = volume_total + volume;
    volume_concave = volume_concave + abs(((1/12)*dot(C-A,cross((D-A),(B-A)))));
    %volume4 = (1/6)*dot(B-R,cross((B-D),(C-A))) - (1/12)*dot(B-A,cross((C-A), (D-A))) + (1/6)*dot(B-A,cross((C-A),(B-D)));
    %volume7 = (1/6)*dot(B-R,cross((B-D),(C-A))) + (1/12)*dot(D-B,cross((A-B), (C-D))) ;
    tetra = [R';A';B';C';D'];
    volume2 = getVol(tetra,[2,3,4;1,3,2;1,4,3;1,2,4]) + getVol(tetra,[2,4,5;1,4,2;1,5,4;1,2,5]);
    %volume3 = getVol(tetra,[2,3,5;1,3,2;1,5,3;1,2,5]) + getVol(tetra,[3,4,5;1,4,3;1,5,4;1,3,5]);
    volume_convex = volume_convex + (1/6)*dot(A-R,cross((B-D),(C-A)));
end
%ratio = (volume_total - volume_convex) / volume_total;
ratio = volume_concave / volume_total;
end
function BuildSubdivide()

planes = [0.0846587867  , 0.0777781233 , 0.993369758   ,2818.81665   ; ... 
-0.0846587867 ,-0.0777781159 ,-0.993369758   ,-2825.77002  ; ...
-0.857492924 , 0.514495790  , -0.00000000    ,1093.73218   ; ...
0.857492924  ,-0.514495790  , 0.00000000     ,-1194.91626  ; ...
0.494009405  , 0.869456530  , -0.00000000    ,-1463.70044 ; ... 
-0.494009405 , -0.869456530 , 0.00000000     ,1370.43140  ];

planes_717 = [0.0846587867  , 0.0777781233 , 0.993369758   ,2818.81665   ; ... 
-0.0846587867 ,-0.0777781159 ,-0.993369758   ,-2825.77002  ; ...
-0.857492924 , 0.514495790  , -0.00000000    ,1093.73218   ; ...
0.857492924  ,-0.514495790  , 0.00000000     ,-1194.91626  ; ...
0.494009405  , 0.869456530  , -0.00000000    ,-1463.70044 ; ... 
-0.494009405 , -0.869456530 , 0.00000000     ,1370.43140  ];

planes = [0.0846587867  , 0.0777781233 , 0.993369758   ,2818.81665   ; ... 
-0.0846587867 ,-0.0777781159 ,-0.993369758   ,-2825.77002  ; ...
-0.857492924 , 0.514495790  , -0.00000000    ,1093.73218   ; ...
0.857492924  ,-0.514495790  , 0.00000000     ,-1194.91626  ; ...
0.494009405  , 0.869456530  , -0.00000000    ,-1463.70044 ; ... 
-0.494009405 , -0.869456530 , 0.00000000     ,1370.43140  ];

pts = [1792.50000 , 665.000000 , -3042.14990 ; ...
1704.50000 , 715.000000 , -3039.18994 ; ...
1744.50000 , 585.000000 , -3032.41992 ; ...
1656.50000 , 635.000000 , -3028.20996 ; ...
1792.50000 , 665.000000 , -3049.14990 ; ...
1704.50000 , 715.000000 , -3046.18994 ; ...
1744.50000 , 585.000000 , -3039.41992 ; ...
1656.50000 , 635.000000 , -3035.20996 ];

pts_717 = [1792.50000 , 665.000000 , -3042.14990 ; ...
1704.50000 , 715.000000 , -3039.18994 ; ...
1744.50000 , 585.000000 , -3032.41992 ; ...
1656.50000 , 635.000000 , -3028.20996 ; ...
1792.50000 , 665.000000 , -3049.14990 ; ...
1704.50000 , 715.000000 , -3046.18994 ; ...
1744.50000 , 585.000000 , -3039.41992 ; ...
1656.50000 , 635.000000 , -3035.20996 ];

pts_61342 = [441.140625 ,-1783.25000 , -3087.29004 ; ...
389.921875 ,-1695.75000 , -3084.12012 ; ...
354.578125 ,-1837.25000 , -3081.08008 ; ...
303.359375 ,-1749.75000 , -3072.65991 ; ...
441.140625 ,-1783.25000 , -3094.29004 ; ...
389.921875 ,-1695.75000 , -3085.14990 ; ...
354.578125 ,-1837.25000 , -3086.21997 ; ...
303.359375 ,-1749.75000 , -3072.65991 ];


planes_61342 = [0.104383603  , -0.00471933465 , 0.994525969    , 3017.22705 ;...
-0.144522741 , 0.0435706489   , -0.988541782  , -2918.47021;...
-0.529282808 , 0.848445415    , -0.00000000 , 1645.13025 ;...
0.529282928  , -0.848445415   , 0.00000000  , -1746.47852;...
0.863017738  , 0.505173564    , -0.00000000  , 520.138550 ;...
-0.863017797 , -0.505173624   , 0.00000000  , -622.122986];
%[k2,dualPoints,volume] = findInterior(planes_717,[0,0,0],false);


planes_645 = [-0.00963814370,  0.0100660194,  0.999902964  , 3006.07764   ;...
-0.0498249233 ,  0.0212220829,  -0.998532474 , -2986.74146  ;... 
-0.882352948  ,  0.470588237 ,  0.00000000   , 268.970581   ;...
0.857492924   ,  -0.514495790,  0.00000000   , -291.976318  ;...
0.494009405   ,  0.869456530 ,  -0.00000000  , -1549.75684  ;...
-0.413802981  ,  -0.910366535,  0.00000000   , 1421.20630 ];

pts_645 = [1051.50000 , 1185.00000 , -3008.67993 ; ...
963.500000 , 1235.00000 , -3008.94995 ; ...
1003.50000 , 1105.00000 , -3007.33008 ; ...
915.500000 , 1145.00000 , -3009.61011 ; ...
1051.50000 , 1185.00000 , -3015.67993 ; ...
963.500000 , 1235.00000 , -3015.94995 ; ...
1003.50000 , 1105.00000 , -3020.33008 ; ...
915.500000 , 1145.00000 , -3009.61011 ];



pts_645_pituba = [2805.50000 225.000000 -3190.96997 ; ...
2717.50000 275.000000 -3181.92993 ; ...
2757.50000 145.000000 -3176.27002 ; ...
2669.50000 185.000000 -3167.50000 ; ...
2805.50000 225.000000 -3197.96997 ; ...
2717.50000 275.000000 -3188.92993 ; ...
2757.50000 145.000000 -3183.27002 ; ...
2669.50000 185.000000 -3174.50000 ];

pts_85140 = [1549.50000 , -2275.00000 , -3094.47998 ; ...
1461.50000 , -2225.00000 , -3091.36011 ;...
1501.50000 , -2365.00000 , -3091.16992 ;...
1413.50000 , -2315.00000 , -3089.90991 ;...
1549.50000 , -2275.00000 , -3094.47998 ;...
1461.50000 , -2225.00000 , -3098.98999 ;...
1501.50000 , -2365.00000 , -3091.16992 ;...
1413.50000 , -2315.00000 , -3099.42993 ];



planes_85140 = [0.0306113940 , 0.0101037854  , 0.999480307  , 3067.96069   ; ...
0.0486132689 , -0.0418346338 , -0.997941136 , -3257.67432  ; ...
-0.882352948 , 0.470588237   , -0.00000000  , 2336.61768   ; ...
0.882352948  , -0.470588237  , 0.00000000   , -2437.79395  ; ...
0.494009405  , 0.869456530   , -0.00000000  , 1212.54614   ; ...
-0.494009376 , -0.869456530  , 0.00000000   , -1314.50964 ];


planes_645_pituba = [0.144395366 , 0.0876506120  , 0.985630453  , 2720.01465   ;...  
-0.144395366, -0.0876506120 , -0.985630453 , -2726.91406  ;...
-0.882352948, 0.470588237   , -0.00000000  , 2268.38232   ;...
0.857492924 , -0.514495790  , 0.00000000   , -2289.93481  ;...
0.494009405 , 0.869456530   , -0.00000000  , -1581.57117  ;...
-0.413802952, -0.910366535  , 0.00000000   , 1273.06482  ];

pts_61155 = [2599.82812 -1804.75000 -3177.69995 ; ...
2511.82812 -1754.75000 -3186.94995 ; ...
2551.82812 -1894.75000 -3176.62012 ; ...
2463.82812 -1844.75000 -3184.40991 ; ...
2599.82812 -1804.75000 -3177.69995 ; ...
2511.82812 -1754.75000 -3188.39990 ; ...
2551.82812 -1894.75000 -3176.62012 ; ...
2463.82812 -1844.75000 -3184.40991 ];



planes_61155 = [-0.0290708691, 0.0799951255 , 0.996371329 , 3377.09082  ; ...  
0.0267062876 , -0.0895394310, -0.995625079, -3388.92065 ; ... 
-0.529282928 , 0.848445415  , -0.00000000 , 2559.22266  ; ... 
0.529282928  , -0.848445415 ,  0.00000000 , -2660.57080 ; ... 
0.863017678  , 0.505173624  , -0.00000000 , -694.591309 ; ... 
-0.863017678 , -0.505173624 , 0.00000000  , 592.606995 ];

planes_61341 = [-0.162805274 0.0885802582  0.982673824  2486.82764  ; ...  
0.162805274  -0.0885802582 -0.982673824 -2493.70654 ; ...
-0.882352948 0.470588237   -0.00000000  -3744.85303 ; ...
0.882352948  -0.470588237  0.00000000   3644.55884  ; ...
0.498283893  0.867013931   -0.00000000  -269.671143 ; ...
-0.498283893 -0.867013931  0.00000000   167.722168  ];

pts_61341 = [-3034.50000 , 2055.00000 , -3220.27002 ; ...
-3121.50000 , 2105.00000 , -3235.96997 ; ...
-3082.50000 , 1965.00000 , -3216.87012 ; ...
-3169.50000 , 2015.00000 , -3239.05005 ; ...
-3034.50000 , 2055.00000 , -3227.27002 ; ...
-3121.50000 , 2105.00000 , -3242.96997 ; ...
-3082.50000 , 1965.00000 , -3223.87012 ; ...
-3169.50000 , 2015.00000 , -3246.05005 ];


pts_68894 = [-3082.50000 1965.00000 -3230.87012; ... 
-3169.50000 2015.00000 -3249.94995; ... 
-3130.50000 1885.00000 -3226.37988; ... 
-3217.50000 1925.00000 -3243.17993; ... 
-3082.50000 1965.00000 -3237.87012; ... 
-3169.50000 2015.00000 -3249.94995; ... 
-3130.50000 1885.00000 -3226.37988; ... 
-3217.50000 1925.00000 -3243.17993];

planes_68894 = [-0.130633935 0.138564438  0.981699944  2496.99561  ; ... 
0.0840663910 -0.154751629 -0.984370232 -2622.52686 ; ...
-0.857492924 0.514495790  -0.00000000  -3749.38770 ; ...
0.857492924  -0.514495790 0.00000000   3654.20605  ; ...
0.498283893  0.867013931  -0.00000000  -167.722168 ; ...
-0.498283893 -0.867013931 0.00000000   74.4436035];

planes_32215 = [-0.0541910902 -0.0727549419  0.995876491 2973.36621     ; ...
0.0541910902  0.0727549419   -0.995876491 -2980.33740   ; ...
-0.857492924  0.514495790    -0.00000000  -3022.23364   ; ...
0.882352948   -0.470588237   0.00000000   2956.02930    ; ...
0.413802952   0.910366535    -0.00000000  328.766479    ; ...
-0.494009405  -0.869456530   0.00000000   -697.689514   ];
pts_32215 = [-2851.50000 935.000000 -3071.89990 ; ...
-2939.50000 975.000000 -3074.97998 ; ...
-2899.50000 845.000000 -3082.38989 ; ...
-2987.50000 895.000000 -3082.25000 ; ...
-2851.50000 935.000000 -3078.89990 ; ...
-2939.50000 975.000000 -3081.97998 ; ...
-2899.50000 845.000000 -3089.38989 ; ...
-2987.50000 895.000000 -3089.25000];

pts_28397 = [-2803.50000 , 1015.00000 , -3060.35010 ; ... 
-2891.50000 , 1065.00000 , -3065.64990 ; ... 
-2851.50000 , 935.000000 , -3064.89990 ; ... 
-2939.50000 , 975.000000 , -3067.97998 ; ... 
-2803.50000 , 1015.00000 , -3067.35010 ; ... 
-2891.50000 , 1065.00000 , -3072.64990 ; ... 
-2851.50000 , 935.000000 , -3071.89990 ; ... 
-2939.50000 , 975.000000 , -3074.97998 ];
planes = CreateHalfs(pts_28397);

planes_28397 = [-0.0529094264 -0.0111526381 0.998537064  2919.43213  ; ... 
0.0529094264  0.0111526381  -0.998537064 -2926.42188 ; ... 
-0.882352948  0.470588237   -0.00000000  -3052.50000 ; ... 
0.857492924   -0.514495790  0.00000000   2926.19458  ; ...
0.494009405   0.869456530   -0.00000000  502.456970  ; ... 
-0.413802952  -0.910366535  0.00000000   -328.766479 ];

pts_degFace = 1.0e+03 * [[0.187500000000000  -0.385000000000000  -3.084879882812500];...
              [0.099500000000000  -0.345000000000000  -3.084870117187500];...
              [0.139500000000000  -0.475000000000000  -3.089310058593750];...
              [0.051500000000000  -0.425000000000000  -3.093459960937500];...
              [0.187500000000000  -0.385000000000000  -3.091290039062500];...
              [0.099500000000000  -0.345000000000000  -3.091870117187500];...
              [0.139500000000000  -0.475000000000000  -3.089310058593750];...
              [0.051500000000000  -0.425000000000000  -3.093459960937500]];
planes_degFace = 1.0e+03 * [[-0.000048481005709  -0.000050230320221   0.000997560277385   3.065920448963012];...
                       [ 0.000019918862483  -0.000012616599804  -0.000999721991519  -3.098144654066119];...
                       [-0.000857492925713   0.000514495755428                   0   0.262821581730895];...
                       [ 0.000882352941176  -0.000470588235294                   0  -0.346617647058824];...
                       [ 0.000413802944301   0.000910366477463                   0   0.272903041766631];...
                       [-0.000413802944301  -0.000910366477463                   0  -0.374698566064722]];

pts_15382 = [[ -344.500000 2185.00000 -3147.56006] ; ...
[ -432.500000 2225.00000 -3158.59009] ; ...
[ -392.500000 2095.00000 -3152.87012] ; ...
[ -480.500000 2145.00000 -3166.83008] ; ...
[ -344.500000 2185.00000 -3147.56006] ; ...
[ -432.500000 2225.00000 -3158.59009] ; ...
[ -392.500000 2095.00000 -3154.64990] ; ...
[ -480.500000 2145.00000 -3166.83008]];

planes = CreateHalfs(pts_degFace);

pts = pts_degFace;
planes = planes;
[k2,dualPoints,vol,bool,center] = findInterior(planes)

figure
hold on
axis off
view(-70,43);
showHalfSpace(pts);
showSurface(pts);


plotPlaness(planes(1:2,:));
bb = getBB(pts);
axis(bb);
view(-153,29);
for i = 1:size(dualPoints,1)
    plot3(dualPoints(i,1),dualPoints(i,2),dualPoints(i,3),'o','MarkerFaceColor','green','MarkerSize',5);
end
new_coords = move(pts,planes);
%exportgraphics(gca,strcat('surface_85140','.png'),'Resolution',333);

%fname = 'final';
%fig2u3d(gca, fname, '-pdf')
%exportgraphics(gca,strcat('cell4','.png'),'Resolution',1000);

%  figure
%  hold on
%  view(-48,71);
%  %showSurface(pts)
%  %showCell(pts)
%  axis off 
%  exportgraphics(gca,strcat('cell2','.png'),'Resolution',1000);
%   hold off
[vol_exact,vol_concave,ratio_init,volume_convex_init] = realVolume(pts);

%[k2,dualPoints,volume_aprox] = findInterior(planes_85140,[0,0,0],false);

out = ratio_vol(pts);

[solid1,solid2,solid3,solid4,normals] = childrenCells(pts);



figure
hold on
axis off
view(110,69);
lam = 5;
%showSurface(pts)
showSurface(solid1 + lam*normals(1,:));
[vol1,vol_concave,ratio_1,vol_convex1] = realVolume(solid1);
ratio_1 = vol_concave/vol1;
plane_s1 = [0.0369095169  0.0119012594  0.999247789  , 3061.92017  ; ...
0.0359498970  -0.0454859063 -0.998317897 , -3248.22559 ; ...
-0.882352948  0.470588237   0.00000000   , 2387.20581  ; ...
0.882352948   -0.470588237  -0.00000000  , -2437.79419  ; ...
0.494009405   0.869456530   0.00000000   , 1212.54602  ; ...
-0.494009435  -0.869456589  0.00000000   , -1263.52795 ];
%plotPlaness(plane_s1(3:4,:));
%bb = getBB(pts);
%axis(bb);
%[k2,dualPoints,volume_aprox] = findInterior(plane_s1,[0,0,0],false);

showSurface(solid2 + lam*normals(2,:) );
[vol2,vol_concave,ratio_2,vol_convex2] = realVolume(solid2);
showSurface(solid3 + lam*normals(3,:));
[vol3,vol_concave,ratio_3,vol_convex3] = realVolume(solid3);
showSurface(solid4 + lam*normals(4,:));
[vol4,vol_concave,ratio_4,vol_convex4] = realVolume(solid4);
%exportgraphics(gca,strcat('cell3','.png'),'Resolution',1000);
%exportgraphics(gca,strcat('partition_2_85140','.png'),'Resolution',333);
volt_tot = vol1+vol2+vol3+vol4;
vol_aprox = vol_convex1  + vol_convex2 + vol_convex3 + vol_convex4;
new_error = (vol_aprox - vol_exact)/vol_exact;
ratio_init
new_error
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
%plot3(R(1),R(2),R(3),'o','MarkerSize',10,'MarkerFaceColor','b');
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
%hold off
end
function new_pts = Subdivide(face)
A = face(1,:)';
B = face(2,:)';
C = face(3,:)';
D = face(4,:)';
new_pts = zeros(5,3);
ksi_midle = [0.5,1,0.5,0,0.5];
eta_midle = [0,0.5,1.0,0.5,0.5];
x = A(1) + ksi_midle.*(B(1) - A(1)) + eta_midle.*(D(1) - A(1)) + ksi_midle.*eta_midle.*( (C(1)-A(1)) - (B(1)-A(1)) - (D(1) - A(1)));
y = A(2) + ksi_midle.*(B(2) - A(2)) + eta_midle.*(D(2) - A(2)) + ksi_midle.*eta_midle.*((C(2)-A(2)) - (B(2)-A(2)) - (D(2) - A(2)));
z = A(3) + ksi_midle.*(B(3) - A(3)) + eta_midle.*(D(3) - A(3)) + ksi_midle.*eta_midle.*((C(3)-A(3)) - (B(3)-A(3)) - (D(3) - A(3)));
hold on
line([x(1),x(3)],[y(1),y(3)],[z(1),z(3)],'LineStyle','--');
line([x(2),x(4)],[y(2),y(4)],[z(2),z(4)],'LineStyle','--');
for i = 1:5
    new_pts(i,:) = [x(i),y(i),z(i)];
   plot3(x(i),y(i),z(i),'o','MarkerFaceColor','b','MarkerSize',5);
end
end
function normals = CreateNormals(face)
A = face(1,:)';
B = face(2,:)';
C = face(3,:)';
D = face(4,:)';
L1 = B - A;
L2 = C - B;
L3 = D - C;
L4 = A - D;
normals = zeros(4,3);
normals(1,:) = -L1 + L4;
normals(2,:) = L1 - L2;
normals(3,:) = -L3 + L2;
normals(4,:) = L3 - L4;
normals(1,:) = normals(1,:)/norm(normals(1,:));
normals(2,:) = normals(2,:)/norm(normals(2,:));
normals(3,:) = normals(3,:)/norm(normals(3,:));
normals(4,:) = normals(4,:)/norm(normals(4,:));
end
function out = ratio_vol(pts)
top = [pts(1,:);pts(2,:);pts(4,:); pts(3,:)];
bottom = [pts(5,:);pts(6,:);pts(8,:); pts(7,:)];
R = sum(pts,1)'/8;
A = top(1,:)';
B = top(2,:)';
C = top(3,:)';
D = top(4,:)';
volume_convexo_top = ((1/6)*dot(A-R,cross((B-D),(C-A))));
volume_concave_top = ((1/12)*dot(C-A,cross((D-A),(B-A))));
A = bottom(1,:)';
B = bottom(2,:)';
C = bottom(3,:)';
D = bottom(4,:)';
volume_convexo_bottom = ((1/6)*dot(A-R,cross((B-D),(C-A))));
volume_concave_bottom = ((1/12)*dot(C-A,cross((D-A),(B-A))));
out = (volume_concave_top + volume_concave_bottom) / ...
    (volume_convexo_top + volume_concave_top + volume_convexo_bottom + volume_concave_bottom);
end
function [solid1,solid2,solid3,solid4,normals] = childrenCells(pts)
top = [pts(1,:);pts(2,:);pts(4,:); pts(3,:)];
bottom = [pts(5,:);pts(6,:);pts(8,:); pts(7,:)];
normals = CreateNormals(top);
new_pts_top = Subdivide(top);
new_pts_bottom = Subdivide(bottom);
solid1 = [[top(1,:) ; new_pts_top(1,:); new_pts_top(5,:);new_pts_top(4,:)];...
    [bottom(1,:) ; new_pts_bottom(1,:); new_pts_bottom(5,:);new_pts_bottom(4,:)]];


solid2 = [[new_pts_top(1,:); top(2,:); new_pts_top(2,:); new_pts_top(5,:)];...
    [new_pts_bottom(1,:); bottom(2,:); new_pts_bottom(2,:); new_pts_bottom(5,:)]];


solid3 = [[new_pts_top(5,:); new_pts_top(2,:); top(3,:);new_pts_top(3,:)];...
    [new_pts_bottom(5,:); new_pts_bottom(2,:);bottom(3,:); new_pts_bottom(3,:)]];

solid4 = [[new_pts_top(4,:); new_pts_top(5,:);new_pts_top(3,:);top(4,:)];...
    [new_pts_bottom(4,:); new_pts_bottom(5,:);new_pts_bottom(3,:);bottom(4,:)]];

ids = [1,2,4,3,5,6,8,7];
solid1 = solid1(ids,:);
solid2 = solid2(ids,:);
solid3 = solid3(ids,:);
solid4 = solid4(ids,:);
end

function new_coords = move(pts,planes)
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
pts = pts(ids_pts,:);
vert_planes = [left;front;right;back];
ids_planes = [4,5,3,6,1,2];
planes = planes(ids_planes,:);
view(62,80);

[o,d] = intersect3D(-planes(6,4)*planes(6,1:3),-planes(5,4)*planes(5,1:3));
z = d/norm(d);
edges_len = zeros(1,4);
for i = 1:4
    p_i = pts(i,:);
    p_j = pts(i+4,:);
    edge = p_j - p_i;
    edges_len(i) = dot(edge,edge);
    if (dot(edge,edge) == 0)
        edge_dir = cross(vert_planes(i,:),vert_planes(mod(i,4)+1,:));
        p_i = p_j + edge_dir;
        p_j = p_j - edge_dir;
        edge = p_j - p_i;
    end
    p_j = rayPlaneIntersection(p_i,edge,planes(5,:));
    p_i = rayPlaneIntersection(p_i,edge,planes(6,:));
    %plot3(p_j(1),p_j(2),p_j(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    %plot3(p_i(1),p_i(2),p_i(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    new_edge = p_j - p_i;
    if dot(edge,new_edge) > 0
        p_j = rayPlaneIntersection(o,d,planes(i,:));
        p_i = rayPlaneIntersection(o,d,planes(mod(i,4)+1,:));
    end
    %h1 = drawPlan(planes(i,1:3)',planes(i,4),'b');
    %h2 = drawPlan(planes(mod(i,4)+1,1:3)',planes(mod(i,4)+1,4),'b');

    %plot3(p_j(1),p_j(2),p_j(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    %plot3(p_i(1),p_i(2),p_i(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    new_coords(i,:) = p_j;
    new_coords(i + 4,:) = p_i;
    %delete([h1,h2]);
end
%check degenerateFaces
for i = 1:4
    ei = edges_len(i);
    j = mod(i,4)+1;
    ej = edges_len(j);
    n = planes(mod(i,4)+1,1:3);
    if (ei == 0 ) && (ej == 0) && abs(dot(n,z)) < 1e-3
        %denerate face
        p_i = rayPlaneIntersection(o,d,planes(i,:));
        p_j = rayPlaneIntersection(o,d,planes(mod(j,4)+1,:));
        new_coords(i,:) = p_i;
        new_coords(i + 4,:) = p_i;
        new_coords(j,:) = p_j;
        new_coords(j + 4,:) = p_j;
    end
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
%axis (getBB(pts));
%set(gca,'XColor','none','YColor','none','ZColor','none')
% view(30,30);
% n = planes(5,1:3);
% d = planes(5,4)-0.2;
% drawPlan(n',d,'b');
% n = planes(6,1:3);
% d = planes(6,4)-0.2;
% drawPlan(n',d,'b');
% view(-108,44)
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
    d = planes(i,4);
    h1 = drawPlan(n',d,'b');
    n = planes(mod(i,4)+1,1:3);
    d = planes(mod(i,4)+1,4);
    h2 = drawPlan(n',d,'b');
    plot3(pt_up_j(1),pt_up_j(2),pt_up_j(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    plot3(pt_down_j(1),pt_down_j(2),pt_down_j(3),'o','MarkerSize',8,'MarkerFaceColor','b');
    quiver3(pt_up_i(1),pt_up_i(2),pt_up_i(3),u_up(1),u_up(2),u_up(3),'g');
    quiver3(pt_down_i(1),pt_down_i(2),pt_down_i(3),u_down(1),u_down(2),u_down(3),'g');
    delete([h1,h2]);
    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i), 'VerticalAlignment','top','HorizontalAlignment','left');
end
%exportgraphics(gca,strcat('proj_points_2','.png'),'Resolution',500);
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
function planes = CreateHalfs(pts)
faces = cell(6,1);


% faces(1) = {[pts(1,:) ; pts(4, :) ; pts(3,:) ; pts(2,:)]};
% faces(2) = {[pts(8,:) ; pts(5, :) ; pts(7,:) ; pts(6,:)]};
% faces(3) = {[pts(2,:) ; pts(8, :) ; pts(4,:) ; pts(6,:)]};
% faces(4) = {[pts(1,:) ; pts(7, :) ; pts(5,:) ; pts(3,:)]};
% faces(5) = {[pts(5,:) ; pts(2, :) ; pts(6,:) ; pts(1,:)]};
% faces(6) = {[pts(4,:) ; pts(7, :) ; pts(3,:) ; pts(8,:)]};

% faces(1) = {[pts(1,:) ; pts(3, :) ; pts(4,:) ; pts(2,:)]};
% faces(2) = {[pts(5,:) ; pts(7, :) ; pts(8,:) ; pts(6,:)]};
% faces(3) = {[pts(2,:) ; pts(4, :) ; pts(8,:) ; pts(6,:)]};
% faces(4) = {[pts(1,:) ; pts(3, :) ; pts(7,:) ; pts(5,:)]};
% faces(5) = {[pts(1,:) ; pts(2, :) ; pts(6,:) ; pts(5,:)]};
% faces(6) = {[pts(3,:) ; pts(4, :) ; pts(8,:) ; pts(7,:)]};


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

function normal = getNormal(pts)
p0 = pts(1,:);
p1 = pts(2,:);
p2 = pts(3,:);
normal = cross(p1-p0,p2-p0);
len = norm(normal);
if len ~= 0
    normal = normal / len;
end
end