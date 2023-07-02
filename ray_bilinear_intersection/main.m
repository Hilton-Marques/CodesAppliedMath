clc;
clear;
close all;

global tol

tol = 1e-6;


%rayOneCell();

draw = Drawer();


model_A = readmatrix("pituba_well_intersection.txt");
%bb = BBIntersection(model_A);
ids_A = model_A(:,1);
model_A = model_A(1:end,2:end);
tested_cells = (readmatrix("bilinear_tested.txt") + 1)';
intersected_cells = [];
ts = [];
all_ts = [];
%[o,o2] = createRayFromCell(model_A,25012);
o = 1.0e+03 *[-1.503541291727291  -0.144718026978615  -3.027352883619063];
o2 =  1.0e+03 * [-1.775541291727291  -0.224718026978615  -3.037852883619063];
%o = [-2962.10010, 1042.30005, -2886.76416]; % Estudo 1
%o2 = [-2033.5000, 568.700012, -3156.89990]; %Estudo 1

% o = [-195,585,-2969.24512];
% o2 = [-195, 585, -3116.24512];
%
% o = [-195,585,-2886.76416];
% o2 = [-195, 585, -2969.24512];
%
% o = [207.250000 , -455.000000 , -2852.36011];
% o2 = [207.250000 , -455.000000 , -3427.19995];
%
% o = [-1731.50000 , 375.000000 ,-2852.36011];
% o2 = [-1731.50000 , 375.000000 , -3427.19995];

d = o2 - o;
draw.drawLine(o,d);
d = d;
plot3(o(1),o(2),o(3),'o','markerfacecolor','blue','color','blue','markersize',5);
%plot3(o2(1),o2(2),o2(3),'o','markerfacecolor','blue','color','blue','markersize',5);

for i = tested_cells
    cell = reshape(model_A(i,:),3,8)';
    if (i == 5267) % (i == 73226) %
        draw.showSurface(cell, draw.m_red);
        draw.setBB(cell);
    end
    [bool,t,valid_ts] = rayCellIntersection(o,d,cell);
    if bool
        intersected_cells(end + 1) = i;
        ts = [ts, t];
        all_ts = [all_ts, valid_ts];
    end
end
[ts, ids] = sort(ts);
intersected_cells = intersected_cells(ids);
%plot
cells = [];
for tsi = all_ts
    if (tsi >= 1)
        break;
    end
    p = o + tsi*d;
    plot3(p(1),p(2),p(3),'o','MarkerSize',5,'color','#00FF00','markerfacecolor','#00FF00');
end
for i = 1:size(ts,2)
    ti = ts(i);
    if ti >= 1
        break;
    end
    cell = reshape(model_A(intersected_cells(i),:),3,8)';
    draw.showSurface(cell, draw.m_red);
    cells = [cells ; cell];
    draw.setBB(cells);
end

draw.setBB(cells);
keyboard;
function [bool,min_t, valid_ts] = rayCellIntersection(o,d,pts)
bool = false;
planes = buildHalfPlanes(pts);
faces = cell(6,1);
faces(1) = {[ pts(1,:) ; pts(2,:) ; pts(4,:) ; pts(3,:)] };
faces(2) = {[ pts(5,:) ; pts(7,:) ; pts(8,:) ; pts(6,:)] };
faces(3) = {[ pts(2,:) ; pts(6,:) ; pts(8,:) ; pts(4,:)] };
faces(4) = {[ pts(1,:) ; pts(3,:) ; pts(7,:) ; pts(5,:)] };
faces(5) = {[ pts(1,:) ; pts(5,:) ; pts(6,:) ; pts(2,:)] };
faces(6) = {[ pts(3,:) ; pts(4,:) ; pts(8,:) ; pts(7,:)] };


ts = realmax * ones(1,8);
%wireframe
%     figure
%     hold on
%     for i = 3:6
%         face = faces{i};
%         for j = 1:4
%             p0 = face(j,:);
%             p1 = face(mod(j,4) + 1 , :);
%             line([p0(1),p1(1)],[p0(2),p1(2)],[p0(3),p1(3)],'linewidth',2)
%         end
%     end
for i = 1:2
    t = rayBilinearIntersection(o,d,faces{i});
    ts(2*i-1:2*i) = t;
end
for i = 3:6
    face = faces{i};
    err_concave = errConcave(face);
    if err_concave ~= 0
        disp('warning');
    end
    t = rayQuadIntersection(o,d,face,planes(i,:));
    ts(i+2) = t;
end

min_t = min(ts);
if min_t < realmax
    bool = true;
end
valid_ts = [];
for t = ts
    if t < realmax
        valid_ts(end+1) = t;
    end
end
valid_ts = unique(valid_ts);
end

function planes = buildHalfPlanes(pts)
faces = cell(6,1);
faces(1) = {[ pts(1,:) ; pts(2,:) ; pts(4,:) ; pts(3,:)] };
faces(2) = {[ pts(5,:) ; pts(7,:) ; pts(8,:) ; pts(6,:)] };
faces(3) = {[ pts(2,:) ; pts(6,:) ; pts(8,:) ; pts(4,:)] };
faces(4) = {[ pts(1,:) ; pts(3,:) ; pts(7,:) ; pts(5,:)] };
faces(5) = {[ pts(1,:) ; pts(5,:) ; pts(6,:) ; pts(2,:)] };
faces(6) = {[ pts(3,:) ; pts(4,:) ; pts(7,:) ; pts(7,:)] };
planes = zeros(6,4);
for i = 1:6
    face_pts = faces{i};
    planes(i,:) = CellIntersection.fittingPlane(face_pts);
end
opp_ids = [2,1,4,3,6,5];
for i = 1:6
    plane = planes(i,:);
    normal = plane(1:3);
    %correct plane
    if (norm(normal) == 0)
        opp_plane = planes(opp_ids(i),:);
        normal = opp_plane(1:3);
        if (norm(normal) == 0)
            % não tem o que fazer
        else
            max_dist = - realmax;
            id_far = -1;
            plane_pts = faces{i};
            pn = - opp_plane(1:3) * opp_plane(4);
            for j = 1:4
                pi = plane_pts(j,:) - pn;
                dist = abs(dot(pi, normal));
                if dist > max_dist
                    id_far = j;
                    max_dist = dist;
                end
            end
            planes(i,:) = [-normal, dot(plane_pts(id_far,:),normal)];
        end
    end
end
end
function [o,p2] = createRayFromCell(model_A,id)
pts = reshape(model_A(id,:),3,8)';
face1 = [ pts(1,:) ; pts(2,:) ; pts(4,:) ; pts(3,:)];
face2 = [ pts(5,:) ; pts(7,:) ; pts(8,:) ; pts(6,:)];
err_concave = errConcave(face1);
rayOneCell(pts);
if err_concave ~= 0
    disp('warning');
end
p1 = pts(1,:);
p2 = pts(4,:);
d = p2 - p1;
o = p1 - (1*d/norm(d)) - 0.5*[0,0,1];
p2 = o + 2*d;
draw = Drawer();
draw.showSurface(pts, draw.m_red);
%plot3(pts(2,1), pts(2,2),pts(2,3),'o','markersize',10);
%plot3(pts(3,1), pts(3,2),pts(3,3),'o','markersize',10);
draw.drawLine(o,d);
plot3(o(1),o(2),o(3),'o','markerfacecolor','blue','color','blue');
draw.setBB(pts);



end
function rayOneCell(pts)
%% inputs
if nargin == 0
pts = [[-1,-1,1.0];
    [1,-1,2.0];
    [-1,1,2.0];
    [1,1,1];
    [-1,-1,-2.5];
    [1,-1,-1];
    [-1,1,-1];
    [1,1,-2.5]];
end
%% Show
faces = [pts(1,:) ; pts(2,:) ; pts(4,:) ; pts(3,:); pts(5,:) ; pts(7,:) ; pts(8,:) ; pts(6,:); ...
         pts(2,:) ; pts(6,:) ; pts(8,:) ; pts(4,:);  pts(1,:) ; pts(3,:) ; pts(7,:) ; pts(5,:); ...
         pts(1,:) ; pts(5,:) ; pts(6,:) ; pts(2,:);  pts(3,:) ; pts(4,:) ; pts(8,:) ; pts(7,:)];
draw = Drawer();
draw.showSurface(pts, draw.m_red);
%plot3(pts(2,1), pts(2,2),pts(2,3),'o','markersize',10);
%plot3(pts(3,1), pts(3,2),pts(3,3),'o','markersize',10);
p1 = pts(1,:);
p2 = pts(4,:);
d = p2 - p1;
o = p1 - d - 0.5*[0,0,1];
p2 = o + 3*d;
draw.drawLine(o,d);
plot3(o(1),o(2),o(3),'o','markerfacecolor','blue','color','blue');
draw.setBB(pts);
d = d/norm(d);
[bool,t, valid_ts] = rayCellIntersection(o,d,pts);
if t < realmax
    for ts = valid_ts
        p = o + ts*d;
        plot3(p(1),p(2),p(3),'o','MarkerSize',7.5,'color','#00FF00','markerfacecolor','#00FF00');
    end
end

end

function rayBilinearIntersectionDebugger()
%% inputs
o = [0,0,2];
d = -[0,0.5,1];
o = o - d;

pts = [[-1,-1,1.0];
    [1,-1,2.0];
    [-1,1,2.0];
    [1,1,1];
    [-1,-1,-2.5];
    [1,-1,-1];
    [-1,1,-1];
    [1,1,-2.5]];
%% Show
draw = Drawer();
draw.showSurface(pts, draw.m_red);
plot3(pts(2,1), pts(2,2),pts(2,3),'o','markersize',10);
plot3(pts(3,1), pts(3,2),pts(3,3),'o','markersize',10);
draw.drawLine(o,d);
plot3(o(1),o(2),o(3),'o','markerfacecolor','blue','color','blue');
draw.setBB(pts);
%% Solver
A = [pts(1,:); pts(2,:); pts(4,:); pts(3,:)];
B = [pts(5,:); pts(7,:); pts(8,:); pts(6,:)];

q00 = A(1,:);
q01 = A(4,:);
q11 = A(3,:);
q10 = A(2,:);

e11 = q11 - q10;
e00 = q01 - q00;

qn = cross((q10 - q00), (q01 - q11));
q00_o = q00 - o;
q10_o = q10 - o;
a = dot(cross(q00 - o, d), e00);
c = dot(qn , d);
b = dot(cross(q10 - o, d), e11);

b = b - (a + c);
det = b^2 - 4*a*c;
if det < 0
    disp("no solution");
end

det = sqrt(det);
if (c == 0)
    u1 = -a/b;
    u2 = -1;
else
    u1 = (-b - det)/2/c;
    u2 = (-b + det)/2/c;
end
t = realmax;
if (0 <= u1 && u1 <= 1)
    pba = e00 + (e11 - e00)*u1; %pb - pa
    opa = q00_o + (q10_o - q00_o)*u1; % o - pa
    n = cross(d, pba);
    det = dot(n,n);
    n = cross(n, opa);
    t1 = dot(n, pba);
    v1 = dot(n, d);
    if (t1 > 0 && 0 <= v1 && v1 <= det)
        t = t1/det;
        u = u1;
        v = v1/det;
    end
end
if (0 <= u2 && u2 <= 1)
    pba = e00 + (e11 - e00)*u2; %pb - pa
    opa = q00_o + (q10_o - q00_o)*u2; % o - pa
    n = cross(d, pba);
    det = dot(n,n);
    n = cross(n, opa);
    t2 = dot(n, pba);
    v2 = dot(n, d);
    if (t2 > 0 && 0 <= v2 && v2 <= det && t2 < t)
        t = t2/det;
        u = u2;
        v = v2/det;
    end
end

q00 = A(1,:);
q01 = A(2,:);
q11 = A(3,:);
q10 = A(4,:);
x = q00(1) + u*(q01(1) - q00(1)) + v*(q10(1) - q00(1)) + u*v*((q11(1)-q00(1)) - (q01(1)-q00(1)) - (q10(1) - q00(1)));
y = q00(2) + u*(q01(2) - q00(2)) + v*(q10(2) - q00(2)) + u*v*((q11(2)-q00(2)) - (q01(2)-q00(2)) - (q10(2) - q00(2)));
z = q00(3) + u*(q01(3) - q00(3)) + v*(q10(3) - q00(3)) + u*v*((q11(3)-q00(3)) - (q01(3)-q00(3)) - (q10(3) - q00(3)));

plot3(x,y,z,'o','markerfacecolor','green','color','green');
pk = o + t*d;
plot3(pk(1),pk(2),pk(3),'o','markerfacecolor','red','color','red');


end

function t = rayBilinearIntersection(o, d, surface)
t = realmax;
t1 = realmax;
t2 = realmax;
q00 = surface(1,:);
q01 = surface(4,:);
q11 = surface(3,:);
q10 = surface(2,:);

e11 = q11 - q10;
e00 = q01 - q00;

qn = cross((q10 - q00), (q01 - q11));
q00_o = q00 - o;
q10_o = q10 - o;
a = dot(cross(q00 - o, d), e00);
c = dot(qn , d);
b = dot(cross(q10 - o, d), e11);

b = b - (a + c);
det = b^2 - 4*a*c;
if det < 0
    return;
end

det = sqrt(det);
if (c == 0)
    u1 = -a/b;
    u2 = -1;
else
    u1 = (-b - det)/2/c;
    u2 = (-b + det)/2/c;
end
if (0 <= u1 && u1 <= 1)
    pba = e00 + (e11 - e00)*u1; %pb - pa
    opa = q00_o + (q10_o - q00_o)*u1; % o - pa
    n = cross(d, pba);
    det = dot(n,n);
    n = cross(n, opa);
    t = dot(n, pba);
    v1 = dot(n, d);
    if (t > 0 && 0 <= v1 && v1 <= det)
        t1 = t/det;
        u1 = u1;
        v1 = v1/det;
    end
end
if (0 <= u2 && u2 <= 1)
    pba = e00 + (e11 - e00)*u2; %pb - pa
    opa = q00_o + (q10_o - q00_o)*u2; % o - pa
    n = cross(d, pba);
    det = dot(n,n);
    n = cross(n, opa);
    t = dot(n, pba);
    v2 = dot(n, d);
    if (t > 0 && 0 <= v2 && v2 <= det)
        t2 = t/det;
        u2 = u2;
        v2 = v2/det;
    end
end
t = [t1,t2];
end


function volume_concave = errConcave(face)
A = face(1,:);
B = face(2,:);
C = face(3,:);
D = face(4,:);
volume_concave = dot(C-A,cross((D-A),(B-A)));
end
function t = rayQuadIntersection(o,d, pts,n_exact)
global tol;
t = realmax;
normal = zeros(1,3);
for i = 1:4
    p0 = pts(i,:);
    p1 = pts(mod(i,4) + 1 , :);
    p2 = pts(mod(i+1,4) + 1 , :);
    normal_i = cross(p1-p0,p2-p0);
    len = norm(normal_i);
    if len == 0
        continue
    end
    normal_i = normal_i / len;
    normal = normal + normal_i;
end
normal_len = norm(normal);
if normal_len == 0
    return
end
centroide = mean(pts);
%quiver3(centroide(1),centroide(2),centroide(3),normal(1),normal(2),normal(3));

normal = normal/normal_len;
denom = dot(normal,d);
if abs(denom) < tol
    return
end
dd = -dot(normal, pts(1,:));
num = dot(o,normal) + dd;
t = - num / denom;
if t < 0
    t = realmax;
    return;
end
pt = o + t * d;
for i = 1:4
    p0 = pts(i,:);
    p1 = pts(mod(i,4) + 1 , :);
    edge = p1 - p0;
    if (dot(edge,edge) < tol)
        continue;
    end
    v = pt - p0;
    value = dot(normal, cross(edge,v));
    if value < 0
        t = realmax;
        return;
    end
end
end

function index = coordToId(i,j,k)
ni = 83;
nj = 45;
nk = 23;
n = ni*nj*nk;
index = (k-1)*ni*nj + (j-1)*ni + (i-1);
end
