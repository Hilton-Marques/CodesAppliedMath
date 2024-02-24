clc;
clear;
close all;

pts = load("Horizonte3.dat");
%pts = load("Horizonte3.dat");
ptCloud = pointCloud(pts);
gridstep = 0.05;
ptCloudDownSampled = pcdownsample(ptCloud,"gridAverage",gridstep);
depth = 8;
mesh = pc2surfacemesh(ptCloudDownSampled,"poisson",depth);
%surfaceMeshShow(mesh,"VerticesOnly",true);
v = mesh.Vertices;
f = mesh.Faces;
surfaceMeshShow(mesh, "Wireframe",true);
meshes = {struct("vertices",v, "faces",f)};
exportVTK(meshes, "horizonte3_poisson.vtk");

function exportVTK(meshes, filename)
% Open the file for writing
fid = fopen(filename, 'w');

% Write the header
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Unstructured Grid\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% Write the mesh data
num_meshes = numel(meshes);
total_vertices = sum(cellfun(@(x) size(x.vertices, 1), meshes));
total_faces = sum(cellfun(@(x) size(x.faces, 1), meshes));
fprintf(fid, 'POINTS %d double\n', total_vertices);
for i = 1:num_meshes
    vertices = meshes{i}.vertices;
    num_vertices = size(vertices, 1);
    fprintf(fid, '%f %f %f\n', vertices.');
end
fprintf(fid, '\n');
fprintf(fid, 'CELLS %d %d\n', total_faces, total_faces * 4);
offset = 0;
for i = 1:num_meshes
    vertices = meshes{i}.vertices;
    num_vertices = size(vertices, 1);
    faces = meshes{i}.faces;
    num_faces = size(faces, 1);
    cell_data = [repmat(3, num_faces, 1) faces - 1 + offset];
    fprintf(fid, '%d %d %d %d\n', cell_data.');
    offset = offset + num_vertices;
end
fprintf(fid, '\n');
fprintf(fid, 'CELL_TYPES %d\n', total_faces);
for i = 1:num_meshes
    faces = meshes{i}.faces;
    num_faces = size(faces, 1);
    cell_types = repmat(5, num_faces, 1);
    fprintf(fid, '%d\n', cell_types);
end

% Close the file
fclose(fid);
end