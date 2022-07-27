% Used as a standard testing script for PROFILER
clear classes %#ok<CLSCR>
addpath([pwd,'\Code'],'-begin');
load('.\Test Case\Test4star.mat');
g = GalaxyModel(Test4star);
g.simulate(16,0.1);
