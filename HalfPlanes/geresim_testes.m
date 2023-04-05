clc;
clear;
close all;

%Equal_cells
planes = [-0.61013147689880509 0.50544946512442357  0.61013147689880509      -2.5272473256221177  ; ...
-0.57735026918962573 0.57735026918962573  -0.57735026918962573   -0.0000000000000000  ; ...
-0.70710678118654757 0.70710678118654757  0.0000000000000000   -0.0000000000000000  ; ...
0.70710678118654757  -0.70710678118654757  -0.0000000000000000  7.0710678118654755  ; ...
-0.0000000000000000  -1.0000000000000000   -0.0000000000000000   0.0000000000000000  ; ...
0.0000000000000000   1.0000000000000000      0.0000000000000000    -10.000000000000000 ];
cgeom_obj = cgeom();
f = [[0,0,0];...
     [0,0,0];...
     [0,0,10];...
     [0,10,0]];

plane = cgeom_obj.fittingPlane(f);

[k2,dualPoints,volume] = cgeom_obj.findInterior(planes,[0,0,0],true);
