clear classes %#ok<CLSCR>
addpath([pwd,'\Code'],'-begin');
load('.\Test Case\Test2starCollision.mat');
g = GalaxyModel(Test2starCollision);
g.simulate(8,0.1);
