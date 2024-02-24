clc;
clear;
close all;
rng('default'); %for colors
addpath(genpath("./"));
addpath("../my_libs/Geresim_Scene/");
addpath("../my_libs/cgeom/");

[vertices_horizon, faces_horizon] = read_vtk_file("meshes/EMB_COMPLETO_160_145_ts.vtk");
%%To do BFS with weights
%[vertices_horizon, faces_horizon] = read_vtk_file("meshes/result_triangles_extended.vtk");
%find_id_by_coordinate(vertices_horizon, target);
%[vertices_horizon, faces_horizon] = removeDuplicatePoints(vertices_horizon, faces_horizon);
[vertices_fault, faces_fault] = read_vtk_file("meshes/SM_study_FALHA_2_ts_ts.vtk");

%n_vertices = size(vertices_fault,1);

%vertices = [vertices_fault; vertices_horizon];
%faces = [faces_fault; faces_horizon + n_vertices ];

%vtk2obj("meshes/SM_study_FALHA_2_ts_ts.vtk");

horizon = Horizon(vertices_horizon, faces_horizon);
%data = readmatrix("meshes/new_triangles.txt") + 1;
%horizon = Horizon(vertices, faces);
fault = Fault(vertices_fault, faces_fault);

tic
Solver(horizon,fault);

toc
keyboard;

m = ManifoldSurfaceMesh(faces);

G = meshGraph(vertices_horizon,faces);
figure
hold on
axis tight
view(-8,71)
trisurf(faces,vertices_horizon(:,1),vertices_horizon(:,2),vertices_horizon(:,3),'EdgeAlpha',0.6,'FaceColor',[0.2588 0.5216 0.9569],'FaceAlpha',0.7)

patches = ManifPatches(G, vertices_horizon);
%DFS(G);



[points,heds,edges,elements] = buildMesh(vertices_horizon,faces);

function vtk2obj(filename)
  [vertices, faces] = read_vtk_file(filename);

  new_filename = split(filename, '.') ;
  new_filename = strcat(new_filename{1} , '.obj');

  % Open the output file for writing
  fileID = fopen(new_filename, 'w');

  % Write the vertex data to the file
  fprintf(fileID, 'v %f %f %f\n', vertices');

  % Write the face data to the file
  fprintf(fileID, 'f %d %d %d\n', faces');
  fclose(fileID);
end
function [vertices, faces] = read_vtk_file(filename)
% READ_VTK_FILE reads a .vtk file and returns the vertices and faces.
% Inputs:
%   filename - the name of the .vtk file
% Outputs:
%   vertices - a Nx3 matrix of vertex coordinates
%   faces - a Mx3 matrix of indices representing the triangular faces

% Open the file
fid = fopen(filename, 'r');
if (fid == -1)
    error(['Could not open file ', filename]);
end

% Find the header and read the format
header = fgetl(fid);
% if (~strcmp(header, '# vtk DataFile Version 3.0'))
%     error('Invalid .vtk file format');
% end

% format = fgetl(fid);
% if (~strcmp(format, 'ASCII'))
%     error('Only ASCII format is supported');
% end

% Find the keyword 'POINTS'
points = fgetl(fid);
strs = splitchar(points);
str_to_compare = strs{1};
while (~strcmp(str_to_compare, 'POINTS'))
    points = fgetl(fid);
    strs = splitchar(points);
    str_to_compare = strs{1};
end

% Read the number of vertices and allocate space for them
nv = sscanf(points, 'POINTS %d double');
vertices = zeros(nv, 3);

% Read the vertices
for i = 1:nv
    line = fgetl(fid);
    coords = sscanf(line, '%f %f %f');
    vertices(i, :) = coords';
end

% Find the keyword 'POLYGONS'
cells = fgetl(fid);
strs = splitchar(cells);
str_to_compare = strs{1};
while (~strcmp(str_to_compare, 'CELLS'))
    cells = fgetl(fid);
    strs = splitchar(cells);
    str_to_compare = strs{1};
end

% Read the number of faces and allocate space for them
nf = sscanf(cells, 'CELLS %d %d');
ncells = nf(1);
faces = zeros(ncells, 3);

% Read the faces
for i = 1:ncells
    line = fgetl(fid);
    indices = sscanf(line, '%d %d %d');
    faces(i, :) = indices(2:4)' + 1; % VTK indices are 0-based, Matlab indices are 1-based
end

% Close the file
fclose(fid);
end



function [points,heds,edges,elements] = buildMesh(vertices,faces)
 nPts = size(vertices,1);
 nFaces = size(faces,1);
 nEdges = nPts + nFaces - 2;
 nHeds = 2* nEdges; 
 points(nPts) = Point;
 heds = Hed.empty;
 edges = Edge.empty;
 elements(nFaces) = Face;
 for i = 1:nPts
     points(i) = Point(vertices(i,:)', i);    
 end
 for i = 1:nFaces
     nHed = length(heds);
     next = [nHed+2,nHed+3,nHed+1];
     for j = 1:3
         nHed = length(heds);
         inc = [faces(i,j),faces(i,mod(j,3)+1)];
         hed = Hed(inc, nHed+1);
         heds(end+1) = hed;
         pt = points(inc(1));
         pt.setHedStart(hed);         
     end
     elements(i) = Triangle(hed, vertices, coords);
     hedInc = nHed - 1; 
     vertices = points(faces(i,:));
 end
 for pt = points
     for i = 1:length(pt.star)
         pti = pt.star(i);
         hedi = pt.hedStar(i);
         if (hedi.edgeId == -1)
             is_edge_with_two_halfs = false;
             for j = 1:length(pti.star)
                 ptj = pti.star(j);
                 hedj = pti.hedStar(j);
                 if ptj.id == pt.id
                     nEdge = length(edges);
                     edges(end+1) = Edge(hedi,hedj, nEdge + 1);
                     hedi.edgeId = nEdge+1;
                     hedj.edgeId = nEdge+1;
                     is_edge_with_two_halfs = true;
                     break;
                 end
             end
             if (~is_edge_with_two_halfs)
                 nEdge = length(edges);
                 edges(end+1) = Edge(hedi, nEdge + 1);
                 hedi.edgeId = nEdge+1;
             end
         end
     end
 end
end

function G = meshGraph(vertices,faces)
nv = size(vertices,1);
nf = size(faces,1);
e1 = faces(:,1:2);
e2 = faces(:,2:3);
e3 = faces(:,[3,1]);
G = zeros(nv,nv);
for i = 1:nf
    G(e1(i,1),e1(i,2)) = 1;
    G(e1(i,2),e1(i,1)) = 1;

    G(e2(i,1),e2(i,2)) = 1;
    G(e2(i,2),e2(i,1)) = 1;

    G(e3(i,1),e3(i,2)) = 1;
    G(e3(i,2),e3(i,1)) = 1;

end

%spy(G)
G = sparse(G);

end
function strings = splitchar(char)
%remove white spaces from begin
n = size(char,2);
strings = {};
if n == 0
    strings = {' '};
    return;
end
del = ' ';
i = 1;
while (strcmp(char(i), del))
    char(i) = [];
end
n = size(char,2);


curr = '';
for i = 1:n
    if (strcmp(char(i), del))
        strings{end+1} = curr;
        curr = '';
    else
        curr(end +1 ) = char(i);
    end
end
if (~isempty(curr))
    strings{end+1} = curr;
end
end

function [new_vertices, new_faces] = removeDuplicatePoints(vertices, faces)
 [new_vertices, ~, old_points_new_ids] = unique(vertices(:,1:3),'rows');
 new_faces = old_points_new_ids(faces);
end
function id = find_id_by_coordinate(vertices, target)
id = -1;
min = realmax;
for i = 1:size(vertices,1)
    d = norm(vertices(i,:) - target);
    if (d < min)
        id = i;
        min =d ;
    end
end
end

function [vertices, faces] = read_segments_vtk_file(filename)
% READ_VTK_FILE reads a .vtk file and returns the vertices and faces.
% Inputs:
%   filename - the name of the .vtk file
% Outputs:
%   vertices - a Nx3 matrix of vertex coordinates
%   faces - a Mx3 matrix of indices representing the triangular faces

% Open the file
fid = fopen(filename, 'r');
if (fid == -1)
    error(['Could not open file ', filename]);
end

% Find the header and read the format
header = fgetl(fid);
% if (~strcmp(header, '# vtk DataFile Version 3.0'))
%     error('Invalid .vtk file format');
% end

% format = fgetl(fid);
% if (~strcmp(format, 'ASCII'))
%     error('Only ASCII format is supported');
% end

% Find the keyword 'POINTS'
points = fgetl(fid);
strs = splitchar(points);
str_to_compare = strs{1};
while (~strcmp(str_to_compare, 'POINTS'))
    points = fgetl(fid);
    strs = splitchar(points);
    str_to_compare = strs{1};
end

% Read the number of vertices and allocate space for them
nv = sscanf(points, 'POINTS %d double');
vertices = zeros(nv, 3);

% Read the vertices
for i = 1:nv
    line = fgetl(fid);
    coords = sscanf(line, '%f %f %f');
    vertices(i, :) = coords';
end

% Find the keyword 'POLYGONS'
cells = fgetl(fid);
strs = splitchar(cells);
str_to_compare = strs{1};
while (~strcmp(str_to_compare, 'CELLS'))
    cells = fgetl(fid);
    strs = splitchar(cells);
    str_to_compare = strs{1};
end

% Read the number of faces and allocate space for them
nf = sscanf(cells, 'CELLS %d %d');
ncells = nf(1);
faces = zeros(ncells, 2);

% Read the faces
for i = 1:ncells
    line = fgetl(fid);
    indices = sscanf(line, '%d %d %d');
    faces(i, :) = indices(2:indices(1)+1)' + 1; % VTK indices are 0-based, Matlab indices are 1-based
end

% Close the file
fclose(fid);
end
