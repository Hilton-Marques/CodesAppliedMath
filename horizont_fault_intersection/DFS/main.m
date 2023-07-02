close all
clc;
clear;

G = [[0,1,1,0,1,0,0,0];
    [1,0,0,1,0,1,0,0];
    [1,0,0,1,0,0,1,0];
    [0,1,1,0,0,0,0,1];
    [1,0,0,0,0,1,1,0];
    [0,1,0,0,1,0,0,1];
    [0,0,1,0,1,0,0,1];
    [0,0,0,1,0,1,1,0]];
GG = graph(G);
%h = plot(GG,'MarkerSize',8,'NodeColor','blue');
s = [1 2 2 3 3 3 4 5 5 5 8 8];
t = [2 3 4 1 4 5 5 3 6 7 9 10];
G = full(adjacency(digraph(s,t)));
n = size(G,1);
g = cell(n,1);
for u = 1:n
    for j = 1:n
        if (G(u,j) ~= 0)
            g{u,1} = [g{u,1},j];
        end
    end
end
%DFS
DFS(g);