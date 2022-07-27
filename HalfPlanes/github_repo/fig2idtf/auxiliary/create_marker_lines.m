function [vertices, edges, line_colors] = create_marker_lines(h, type)
% File:      create_marker_lines.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.14
% Language:  MATLAB R2012a
% Purpose:   create linesets for markers which are not points
% Copyright: Ioannis Filippidis, 2012-

% depends
%   get_line_xyz, axes_extremal_xyz

p = get_line_xyz(h);

switch type
    case 'o'
        type = 'circle_marker';
    case 'x'
        type = 'x_marker';
    case '+'
        type = 'plus_marker';
    case '*'
        type = 'star_marker';
    case 'square'
        type = 'square_marker';
    case 'diamond'
        type = 'diamond_marker';
    case 'v'
        type = 'down_triangle_marker';
    case '^'
        type = 'up_triangle_marker';
    case '<'
        type = 'left_triangle_marker';
    case '>'
        type = 'right_triangle_marker';
    case 'pentagram'
        type = 'pentagram_marker';
    case 'hexagram'
        type = 'hexagram_marker';
    otherwise
        error('Unknown marker style.')
end
disp(['      MarkerStyle = ', type] );

% line with single point only ?
npoints = size(p, 2);
if npoints == 1
    % then we cannot take the minimum distance between line points....
    r = marker_size_from_bounding_sphere(h);
else
    dp = diff(p, 1, 2);
    d = vnorm(dp);
    r = min(d) /6;  % 1/6 minimal distance between consecutive points,
                    % used as marker size
    
    if r == 0
        msg = ['Minimal distance between consecutive line points = 0.',...
               'Tyring mean...'];
        warning('markers:size', msg)
        
        r = mean(d) /6;
    end
    
    if r == 0
        msg = ['All line points are the same. Mean distance = 0.',...
               'Using figure bounding sphere...'];
        warning('markers:size', msg)
        
        r = marker_size_from_bounding_sphere(h);
    end
end
[marker_vertices, marker_edges] = feval(type, r);
[vertices, edges] = copy_marker_on_points(marker_vertices, marker_edges, p);

linecolor = {get(h, 'Color') };
line_colors = repmat(linecolor, 1, npoints);

function [r] = marker_size_from_bounding_sphere(h)
ax = get(h, 'Parent');
xyz_minmax = axes_extremal_xyz(ax);
v = xyz_minmax(2:2:end) -xyz_minmax(1:2:end);
d = norm(v);
r = d /100; % 1/100 world boounding sphere radius

function [vertices, edges] = copy_marker_on_points(marker_vertices, marker_edges, p)
n = size(p, 2);
m = size(marker_vertices, 2);

marker_vertices = [marker_vertices; zeros(1, m) ];
vertices = cell(1, n);
edges = cell(1, n);
for i=1:n
    curp = p(:, i);
    
    v = bsxfun(@plus, marker_vertices, curp);
    e = marker_edges;
    
    vertices{1, i} = v;
    edges{1, i} = e;
end

function [vertices, edges] = circle_marker(r)
n = 10;
t = linspace(0, 2*pi, n);

vertices = r *[cos(t); sin(t) ];
edges = [1:n; 2:n, 1] -1;

function [vertices, edges] = x_marker(r)
vertices = r *[-1, 1; -1, -1; 1, -1; 1, 1].';
edges = [1, 3; 2, 4].' -1;

function [vertices, edges] = plus_marker(r)
vertices = r *[-1, 0; 1, 0; 0, -1; 0, 1].';
edges = [1, 2; 3, 4].' -1;

function [vertices, edges] = star_marker(r)
vertices = r *[-1, 0; -1, -1; 0, -1; 1, -1; 1, 0; 1, 1; 0, 1; -1, 1].';
edges = [1, 5; 2, 6; 3, 7; 4, 8].' -1;

function [vertices, edges] = square_marker(r)
vertices = r *[-1, 1; -1, -1; 1, -1; 1, 1].';
edges = [1, 2; 2, 3; 3, 4; 4, 1].' -1;

function [vertices, edges] = diamond_marker(r)
vertices = r *[-1, 0; 1, 0; 0, -1; 0, 1].';
edges = [1, 3; 3, 2; 2, 4; 4, 1].' -1;

function [vertices, edges] = down_triangle_marker(r)
vertices = r *[-1, 0; 0, -1; 1, 0].';
edges = [1, 2; 2, 3].' -1;

function [vertices, edges] = up_triangle_marker(r)
vertices = r *[-1, 0; 0, 1; 1, 0].';
edges = [1, 2; 2, 3].' -1;

function [vertices, edges] = left_triangle_marker(r)
vertices = r *[0, -1; 1, 0; 0, 1].';
edges = [1, 2; 2, 3].' -1;

function [vertices, edges] = right_triangle_marker(r)
vertices = r *[0, -1; -1, 0; 0, 1].';
edges = [1, 2; 2, 3].' -1;

function [vertices, edges] = pentagram_marker(r)
r = r /100;
vertices = r *[0, 85; 75, 75; 100, 10; 125, 75;
               200, 85; 150, 125; 160, 190; 100, 150;
               40, 190; 50, 125; 0, 85].';
edges = [1:10, 11; 2:11, 1] -1;

function [vertices, edges] = hexagram_marker(r)
r1 = r;
r2 = r1 /2;

n = 7;
t = linspace(0, 2*pi, n);
x1 = r1 *[cos(t); sin(t) ];

t2 = t +pi /n;
x2 = r2 *[cos(t2); sin(t2) ];

m = 2 *n;
x(:, 1:2:m) = x1;
x(:, 2:2:m) = x2;

vertices = x;
edges = [1:m; 2:m, 1] -1;
function [xyz_minmax, dim] = axes_extremal_xyz(ax)
%AXES_EXTREMAL_XYZ  Minimum and maximum coordinates of graphics objects.
%
% usage
%   [xyz_minmax, dim] = AXES_EXTREMAL_XYZ(ax)
%
% input
%   ax = axes object handle
%
% output
%   xyz_minmax = coordinate extremals of objects in axes ax
%              = [xmin, xmax, ymin, ymax, zmin, zmax]
%   dim = dimension of space, depending on plotted objects
%       = 2 || 3
%
% See also PLOT_SCALINGS.
%
% File:      axes_extremal_xyz.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2009.02.12 - 2012.06.11
% Language:  MATLAB R2012a
% Purpose:   min & max coordinates of objects in plot
% Copyright: Ioannis Filippidis, 2009-

% todo
%   properly process quivergroup and cotourgroup objects

%% input
if nargin < 1
    ax = gca;
elseif isempty(ax) || ~ishandle(ax)
    ax = gca;
end

%% get plot objects, find coordinate limits
obj = get(ax, 'Children');
obj2 = findall(ax, 'Type', 'line'); % get hidden lines

obj = [obj; obj2];

% no graphics ?
if isempty(obj)
    warning('plot_scalings:extremalcoor', 'Axes do not have children.')
    xyz_minmax = [];
    dim = [];
    return
end

%% clear text objects
type = get(obj, 'type');
idx = strcmp(type, 'text');
txtobj = obj(idx);
obj(idx) = [];

%% clear light objects
type = get(obj, 'type');
idx = strcmp(type, 'light');
obj(idx) = [];

[coor_min, coor_max] = nontxt_obj_coor_extremals(obj);

% handle text objects separately
pos = get(txtobj, 'Position');
if iscell(pos)
    pos = cell2mat(pos);
end
mintxtcoor = min(pos, [], 1);
maxtxtcoor = max(pos, [], 1);

% no plot contents ?
if isempty(coor_min) && isempty(mintxtcoor)
    xyz_minmax = [];
    dim = [];
    return
end

if isempty(mintxtcoor)
    xyz_minmax(1, 1:2:6) = coor_min;
    xyz_minmax(1, 2:2:6) = coor_max;
elseif isempty(coor_min)
    xyz_minmax(1, 1:2:6) = mintxtcoor;
    xyz_minmax(1, 2:2:6) = maxtxtcoor;
else
    xyz_minmax(1, 1:2:6) = min(coor_min, mintxtcoor);
    xyz_minmax(1, 2:2:6) = max(coor_max, maxtxtcoor);
end

% dimension ?
if xyz_minmax(1, 5) == xyz_minmax(1, 6)
    xyz_minmax = xyz_minmax(1, 1:4);
end

dim = size(xyz_minmax, 2) /2;
function [coor_min, coor_max] = nontxt_obj_coor_extremals(obj)
% no non-text objects
if isempty(obj)
    coor_max = [];
    coor_min = [];
    return
end

X = get(obj, 'XData');
[minX, maxX] = coor_extremals_wrapper(X);

Y = get(obj, 'YData');
[minY, maxY] = coor_extremals_wrapper(Y);

Z = get(obj, 'ZData');
[minZ, maxZ] = coor_extremals_wrapper(Z);

% ensure [1 x 3] matrices are returned
if isempty(minZ)
    minZ = 0;
    maxZ = 0;
end

coor_min = [minX, minY, minZ];
coor_max = [maxX, maxY, maxZ];

function [minC, maxC] = coor_extremals_wrapper(C)
[minC, maxC] = coor_extremals(C);
function [minC, maxC] = coor_extremals(C)
% usage
%   [minC, maxC] = coor_extremals(C)
%
% input
%   C = cell matrix of nested cell/numeric matrices, with numeric or empty
%       matrices at the bottom nesting level, of coordinates. These belong
%       to different graphics objects.
%
% output
%   minC = minimum coordinate over all graphics objects
%   maxC = maximum coordinate over all graphics objects
%
% See also PLOT_SCALINGS, MIN_CELL, MAX_CELL.

% depends
%   min_cell, max_cell

maxC = max_cell(C);
minC = min_cell(C);
function answer = max_cell(a)
%MAX_CELL   maximum of 2D cell or 2D numeric matrix.
%   MAX_CELL(A) calculates the maximum of the contents of A,
%   provided A is a 2 dimesional numeric matrix or
%   a 2 dimensional cell matrix containing a mix of 
%   2D cell matrices and 2D numeric matrices, which in
%   turn may recursively contain others.
%
%   They should not contain any characters or strings,
%   that is, recursive parsing of the elements should
%   end with 2D numeric matrices, empty numeric matrices
%   or empty cell matrices.
%
%   If A is the empty array [] MAX_CELL returns an empty array.
%
%   See also EXTREMA, MIN_CELL, MAX, ISCELL, ISNUMERIC, ISEMPTY.
%
% File:      max_cell.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2009.07.15 - 2010.02.21
% Language:  MATLAB R2011b
% Purpose:   Maximum of 2D Cell/Numeric Matrix
% Copyright: Ioannis Filippidis, 2009-

if iscell(a)
    if isempty(a)
        answer = [];
    else
        a_maxs = zeros(1,numel(a));
        k = 1;
        for i=1:size(a,1)
            for j=1:size(a,2)
                b = a{i, j};
                temp = max_cell(b);
                if (isempty(temp) == 0)
                    a_maxs(k) = temp;
                    k = k + 1;
                end
            end
        end
        answer = max(a_maxs(1:k-1));
    end
elseif isnumeric(a)
    if isempty(a)
        answer = [];
    else
        n = ndims(a);
        for i=1:n
            a = max(a, [], n+1-i);
        end
        answer = a;
    end
else
    disp('Cannot find max. Neither cell nor numeric matrix!')
end
function answer = min_cell(a)
%MIN_CELL   minimum of 2D cell or 2D numeric matrix.
%   MIN_CELL(A) calculates the minimum of the contents of A,
%   provided A is a 2 dimesional numeric matrix or
%   a 2 dimensional cell matrix containing a mix of 
%   2D cell matrices and 2D numeric matrices, which in
%   turn may recursively contain others.
%
%   They should not contain any characters or strings,
%   that is, recursive parsing of the elements should
%   end with 2D numeric matrices, empty numeric matrices
%   or empty cell matrices.
%
%   If A is the empty array [] MIN_CELL returns an empty array.
%
%   See also EXTREMA, MAX_CELL, MIN, ISCELL, ISNUMERIC, ISEMPTY.
%
% File:      min_cell.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2009.07.15 - 2010.02.21
% Language:  MATLAB R2011b
% Purpose:   Minimum of 2D Cell/Numeric Matrix
% Copyright: Ioannis Filippidis, 2009-

if iscell(a)
    if isempty(a)
        answer = [];
    else
        a_mins = zeros(1,numel(a));
        k = 1;
        for i=1:size(a,1)
            for j=1:size(a,2)
                b = a{i, j};
                temp = min_cell(b);
                if (isempty(temp) == 0)
                    a_mins(k) = temp;
                    k = k + 1;
                end
            end
        end
        answer = min(a_mins(1:k-1));
    end
elseif isnumeric(a)
    if isempty(a)
        answer = [];
    else
        n = ndims(a);
        for i=1:n
            a = min(a, [], n+1-i);
        end
        answer = a;
    end
else
    disp('Cannot find min. Neither cell nor numeric matrix!')
end
