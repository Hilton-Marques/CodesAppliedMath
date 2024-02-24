% A simple file to converter a .vtk to .obj file %
fvtk2obj("result_triangles_0408.vtk");

function fvtk2obj(filename)
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
function obj = readObj(fname)
%
% obj = readObj(fname)
%
% This function parses wavefront object data
% It reads the mesh vertices, texture coordinates, normal coordinates
% and face definitions(grouped by number of vertices) in a .obj file 
% 
%
% INPUT: fname - wavefront object file full path
%
% OUTPUT: obj.v - mesh vertices
%       : obj.vt - texture coordinates
%       : obj.vn - normal coordinates
%       : obj.f - face definition assuming faces are made of of 3 vertices
%
% Bernard Abayowa, Tec^Edge
% 11/8/07
% set up field types
v = []; vt = []; vn = []; f.v = []; f.vt = []; f.vn = [];
fid = fopen(fname);
% parse .obj file 
while 1    
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end  % exit at end of file 
     ln = sscanf(tline,'%s',1); % line type 
     %disp(ln)
    switch ln
        case 'v'   % mesh vertexs
            v = [v; sscanf(tline(2:end),'%f')'];
        case 'vt'  % texture coordinate
            vt = [vt; sscanf(tline(3:end),'%f')'];
        case 'vn'  % normal coordinate
            vn = [vn; sscanf(tline(3:end),'%f')'];
        case 'f'   % face definition
            fv = []; fvt = []; fvn = [];
            str = textscan(tline(2:end),'%s'); str = str{1};
       
           nf = length(findstr(str{1},'/')); % number of fields with this face vertices
           [tok str] = strtok(str,'//');     % vertex only
            for k = 1:length(tok) fv = [fv str2num(tok{k})]; end
           
            if (nf > 0) 
            [tok str] = strtok(str,'//');   % add texture coordinates
                for k = 1:length(tok) fvt = [fvt str2num(tok{k})]; end
            end
            if (nf > 1) 
            [tok str] = strtok(str,'//');   % add normal coordinates
                for k = 1:length(tok) fvn = [fvn str2num(tok{k})]; end
            end
             f.v = [f.v; fv]; f.vt = [f.vt; fvt]; f.vn = [f.vn; fvn];
    end
end
fclose(fid);
% set up matlab object 
obj.v = v; obj.vt = vt; obj.vn = vn; obj.f = f;
end