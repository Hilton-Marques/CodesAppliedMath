clc
clear all;
close all;

addpath('tools');


path = '../data/intel.dat';
%path = '../data/intel.dat';
% 
%poses,landmarks,edges,p,l] = ImportData(path);

%graph = Graph(poses,edges);
%graph.Show();
%G.ErrorDeformation();

% load the graph into the variable g
% only leave one line uncommented

% simulation datasets
%load ../data/simulation-pose-pose.dat
g = importdata(path);
[p, l] = get_poses_landmarks(g);
NLAND = size(l,1);
NPOSE = size(p,1);
DOF = 3;
DIM = 2;
%load ../data/simulation-pose-landmark.dat

% real-world datasets
%load ../data/intel.dat
%load ../data/dlr.dat

% plot the initial state of the graph
%plot_graph(g, 0);

disp( compute_global_error(g));

% the number of iterations
numIterations = 100;

% maximum allowed dx
EPSILON = 10^-4;

% Error
err = 0;

% carry out the iterations
e = [];
global flag;
flag = true;
for i = 1:numIterations

  [dx,error] = linearize_and_solve(g);
  e(end+1) = error;
  % TODO: apply the solution to the state vector g.x
  if flag
      for k = 1:NLAND
          id = l(k);
          g.x(id:id+DIM-1,1) = g.x(id:id+DIM-1,1) + dx(id:id+DIM-1,1);
      end      
      for k = 1:NPOSE
          id = p(k);
          g.x(id:id+DOF-1,1) = t2v(exp(dx(id:id+DOF-1,1),v2t(g.x(id:id+DOF-1,1))));
      end
  else
      
      g.x = g.x + dx;
  end
  
  [p, l] = get_poses_landmarks(g);
  poses = [];
  landmarks = [];
  if (~isempty(l))
      landmarkIdxX = l+1;
      landmarkIdxY = l+2;
      landmarks = [g.x(landmarkIdxX), g.x(landmarkIdxY)];
  end
  if (~isempty(p) > 0)
      pIdxX = p+1;
      pIdxY = p+2;
      poses = [g.x(pIdxX), g.x(pIdxY)];
  end
  %graph.UpdatePose(poses);
  %graph.Show();
  % plot the current state of the graph
  %plot_graph(g, i);

  err = compute_global_error(g);

  % Print current error
  disp(err);

  % TODO: implement termination criterion as suggested on the sheet
  if (max(abs(dx)) < EPSILON)
    break;
  end
end
plot_graph(g, i)
plot(e);
keyboard;
disp(i)
%graph.Save();
function Q = exp(v,P)
theta = v(3);
t = v(1:2);

%Rotation part
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

%Translation part
V = eye(2);
if (theta ~= 0)
    V = [sin(theta), -(1-cos(theta)); 1-cos(theta), sin(theta)] .* (1/theta);
end
t = V * t;

Q = P*[R, t; 0 0 1];
end
function [poses,landmarks,edges,p,l] = ImportData(path)
g = importdata(path);
[p, l] = get_poses_landmarks(g);
poses = [];
landmarks = [];
if (~isempty(l))
    landmarkIdxX = l+1;
    landmarkIdxY = l+2;
    landmarks = [g.x(landmarkIdxX), g.x(landmarkIdxY)];
end
if (~isempty(p) > 0)
    pIdxX = p+1;
    pIdxY = p+2;
    poses = [g.x(pIdxX), g.x(pIdxY)];
end
ids = [pIdxX,pIdxY];
poseEdgesP1 = [];
poseEdgesP2 = [];
landmarkEdgesP1 = [];
landmarkEdgesP2 = [];
ids_from = [];
ids_to = [];
for eid = 1:length(g.edges)
    edge = g.edges(eid);
    if (strcmp(edge.type, 'P') ~= 0)
        poseEdgesP1 = [poseEdgesP1, g.x(edge.fromIdx:edge.fromIdx+1)];
        poseEdgesP2 = [poseEdgesP2, g.x(edge.toIdx:edge.toIdx+1)];
        ids_from = [ids_from;[edge.fromIdx:edge.fromIdx+1']];
        ids_to = [ids_to;[edge.toIdx:edge.toIdx+1']];
    elseif (strcmp(edge.type, 'L') ~= 0)
        landmarkEdgesP1 = [landmarkEdgesP1, g.x(edge.fromIdx:edge.fromIdx+1)];
        landmarkEdgesP2 = [landmarkEdgesP2, g.x(edge.toIdx:edge.toIdx+1)];
    end
end
id_tot = 1:length(ids);
[ids_i,ids_ii] = ismember(ids_from,ids,'rows');
[ids_j,ids_jj] = ismember(ids_to,ids,'rows');
edges = [id_tot(ids_ii)',id_tot(ids_jj)'];

end
