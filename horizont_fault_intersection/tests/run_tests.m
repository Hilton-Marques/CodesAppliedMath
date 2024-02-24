clc;
clear;
close all;

addpath(genpath("../"));
addpath("../..//my_libs/Geresim_Scene/");
addpath("../../../my_libs/cgeom/");
addpath(genpath("../../../gptoolbox/"));


filename = "../meshes/EMB_COMPLETO_160_145_ts.vtk";
filename = "neigh_horizon_complicated.vtk";
tests = Tests(filename);
boundary_id  = 1;
tests.TestExtension(boundary_id);

%------ Mesh simplification ----- %
filename = "cube.obj";
tests.MeshDecimation(filename);

%------- Normal Boundary --------%
vertex_id = 442; %with problem
vertex_id = 52;
tests.TestNormal(vertex_id);

tests.TestNormal(vertex_id);
vertex_id = 1087;
tests.TestBoundaryNormal(vertex_id);

%-------- Remove Triangles of small area ----------%
eps = 0.000015;
FF_2 = tests.RemoveBadTriangle(eps);
FF_1 = tests.MyRemoveBadTriangle(eps);

%------ SDF offset -------------%
tests.TestRetriangulation();
boundary_id  = 2;
tests.TestSDFExtension(boundary_id);

boundary_id  = 2;
tests.TestExtension(boundary_id);





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
