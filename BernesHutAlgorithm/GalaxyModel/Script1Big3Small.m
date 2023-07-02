clear classes %#ok<CLSCR>
addpath([pwd,'\Code'],'-begin');
load('.\Test Case\Test1Big3Small.mat');
g = GalaxyModel(Test1Big3Small);
g.simulate(40,0.1);
