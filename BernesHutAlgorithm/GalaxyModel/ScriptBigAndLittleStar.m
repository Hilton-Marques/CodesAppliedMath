clear classes %#ok<CLSCR>
addpath([pwd,'\Code'],'-begin');
load('.\Test Case\TestBigAndLittleStar.mat');
g = GalaxyModel(TestBigAndLittleStar);
g.simulate(40,0.05);
