clc;
clear;
close all;

addpath(genpath("./"));

[vertices, faces] = read_vtk_file("meshes/EMB_COMPLETO_160_145_ts.vtk");
[points,heds,edges,elements] = buildMesh(vertices,faces);

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
strs = strsplit(points);
str_to_compare = strs{1};
while (~strcmp(str_to_compare, 'POINTS'))
    points = fgetl(fid);
    strs = strsplit(points);
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
strs = strsplit(cells);
str_to_compare = strs{1};
while (~strcmp(str_to_compare, 'CELLS'))
    cells = fgetl(fid);
    strs = strsplit(cells);
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
         hedj = Hed([faces(i,j),faces(i,mod(j,3)+1)],nHed+1,i,next(j));
         heds(end+1) = hedj;
         pt = points(faces(i,j));
         ptStar = points(faces(i,mod(j,3)+1));
         pt.star(end+1) = ptStar;
         pt.hedStar(end+1) = hedj;
     end
     hedInc = nHed - 1; 
     vertices = points(faces(i,:));
     elements(i) = Triangle(hedInc, heds, vertices);        
 end
 for pt = points
     for i = 1:length(pt.star)
         pti = pt.star(i);
         hedi = pt.hedStar(i);
         if (hedi.edgeId == -1)
             for j = 1:length(pti.star)
                 ptj = pti.star(j);
                 hedj = pti.hedStar(j);
                 if ptj.id == pt.id
                     nEdge = length(edges);
                     edges(end+1) = Edge(hedi,hedj, nEdge + 1);
                     hedi.edgeId = nEdge+1;
                     hedj.edgeId = nEdge+1;
                 end
             end
         end
     end
 end
end
