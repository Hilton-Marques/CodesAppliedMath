clear classes %#ok<CLSCR>
addpath([pwd,'\Code'],'-begin');
load('.\Test Case\Test2star.mat');
g = GalaxyModel(Test2star);
g.simulate(20,0.1);
