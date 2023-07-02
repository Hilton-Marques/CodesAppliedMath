clear classes %#ok<CLSCR>
addpath([pwd,'\Code'],'-begin');
load('.\Test Case\TestBigMiddleSmall.mat');
g = GalaxyModel(TestBigMiddleSmall);
g.simulate(630,0.5);

