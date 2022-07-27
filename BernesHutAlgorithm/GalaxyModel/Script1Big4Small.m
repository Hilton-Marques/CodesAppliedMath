clear classes %#ok<CLSCR>
addpath([pwd,'\Code'],'-begin');
load('.\Test Case\Test1Big4Small.mat');
g = GalaxyModel(Test1Big4Small);
g.simulate(40,0.1);
