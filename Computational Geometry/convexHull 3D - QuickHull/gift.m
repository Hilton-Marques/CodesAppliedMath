clc;
clear;
close all;

points = [1.0 , 0.0 , 0.0;... 
          2.0 , 0.0 , 0.0;... 
          1.0 , 1.0 , 0.0;...
          1.0 , 0.0 , 1.0];
      

% 
% points = [0.164082691  ,-0.0418397672 ,0.855258405   ;...
% -0.0329033844 ,0.00839009788 , -0.171504349;...
% -0.769165754  ,	0.410221756	,	-0.00000000 ;...
% 0.00908060931 ,-0.00544836558, 0.00000000 ;...
% 0.00529624289 ,0.00921546202,	-0.00000000 ;...
% -0.0538823903 ,-0.117194198 ,	0.00000000 ;...
% 0.00242908997 ,	0.00129536062,  0.170865238 ;...
% -0.0123904785 ,-0.00660746964, -0.87156176;...
% -0.00535226800,0.00857973564,-0.00000000;...
% 0.461387157 	,-0.739607930 ,	0.00000000 ;...
% 0.0267334040 ,	0.0155590735 ,	-0.00000000;...
% -0.0124093071 ,-0.00722232461, 0.00000000 ];

points = [-0.0075604319745329078, -0.0020193218133660981, 0.1414946076038145;...
1.0015830866495388, 0.91305205653544619, -11.646957044215556;...
-0.73171360880526748, 0.38860918628704361, -0.0000000000000000;...
0.0090507332936479515, -0.0054089909123961836, 0.0000000000000000;...
0.0044512432529275094, 0.0098382557244373978, -0.0000000000000000;...
-4.8334374573500867, -10.682991477091276, 0.0000000000000000;...
0.0015143759955473066, 0.022788195872557430, 0.15348288214910985;...
-0.020357382027326384, -0.30633608196063233, -2.0632324308834025;...
-10.054572649157180, 6.0327440088325073, -0.0000000000000000;...
0.0091540606631215528, -0.0048821656663856792, 0.0000000000000000;...
4.8520655085148210, 10.674544607960506, -0.0000000000000000;...
-0.0053010469675394342, -0.0093298424070305339, 0.0000000000000000];
figure
hold on
view(30,30);

dots = zeros(size(points,1));
for i = 1:size(points,1)
    pi = points(i,:);
    normi = norm(pi);
    %quiver3(0,0,0,pi(1),pi(2),pi(3));
    for j = 1:size(points,1)
        pj = points(j,:);
        normj = norm(pj);
        dots(i,j) = norm(pi - pj);
    end
end
hold off;

%points = xlsread('input','Planilha1');
%points = makeRand(1000);

figure
hold on
view(30,30);
tic
[k1,av1] = convhull(points);
toc
trisurf(k1,points(:,1),points(:,2),points(:,3),'FaceColor','cyan');
axis equal;
hold off

vertices(length(points)) = Verts();

figure
hold on
view(30,30);
for i = 1:length(points)
    p = points(i,:);
    vertices(i) = Verts(i,p);
    plot3(p(1),p(2),p(3),'o');
    text(p(1),p(2),p(3),num2str(i));
end

edges = Edge.empty;
supportEdges = Edge.empty;
triangles = Triangle.empty;
triangle = findTriangleOnHull(vertices);
triangle.show();
supportEdges(end+1) = triangle.edges(1).getTwin();

for i = 2:3
    edgei = triangle.edges(i).getTwin();
    edges(end+1) = edgei;
    supportEdges(end+1) = edgei;
end

triangles(end+1) = triangle;
nel = 1;
count = 0;
while ( length(edges) ~= 0)
    edgei = edges(end);
    edges(end) = [];
    if (~edgei.isProcessed())
        if (nel == 6)
            a = 1;
        end
        q = pivotOnEdge(vertices,edgei);
        
        t = Triangle(edgei,q);
        
        triangles(end+1) = t;
        nel = nel + 1;
        t.id = nel;
        
        for i = 2:3
            flag = false;
            edgej = t.edges(i).getTwin;
            for k = 1:length(supportEdges)
                edgek = supportEdges(k);
                if ((edgej.p0.id == edgek.p0.id) && (edgej.p1.id == edgek.p1.id) || ...
                        (edgej.p1.id == edgek.p0.id) && (edgej.p0.id == edgek.p1.id)  )
                    edgek.processed = true;
                    flag = true;
                    t.edges(i) = edgek;
                    break;
                end
            end
            if ~flag
                edges(end+1) = edgej;
                supportEdges(end+1) = edgej;
            end
        end
        edgei.show();
        t.show();
        %edgei.show();
        count = count + 1;
        edgei.processed = true;
    end
    m = showIds(triangles);
end
toc
figure
hold on
view(30,30)
for t = triangles
    t.show();
end
B = sort(m,2);
[~,id] = sort(B(:,1));
B = B(id,:);
C = sort(k1,2);
[~,id] = sort(C(:,1));
C = C(id,:);
[C,ia,ic] = unique(B,'rows');
C = unique(B,'rows');
keyboard


function points = makeRand(n)
points = zeros(n,3);
for i = 1:n
    points(i,:) = [rand,rand,rand];
end
end

function triangle = findTriangleOnHull(vertices)
edge = findEdgeOnHull(vertices);
p = pivotOnEdge(vertices,edge);
triangle = Triangle(edge, p);

end
function edge = findEdgeOnHull(vertices)
pRight = vertices(1).coord;
x = pRight(1);
id = 1;
for i =  2:length(vertices)
    pi = vertices(i).coord;
    x_ = pi(1);
    if (x_ > x)
        x = x_;
        pRight = pi;
        id = vertices(i).id;
    end
end
suportVertex = Verts(length(vertices) + 1, pRight + [0,1,0]);
edge = Edge ( suportVertex , vertices(id));
r = pivotOnEdge(vertices,edge);
edge = Edge(vertices(id),r);
end
function v = pivotOnEdge(vertices,edge)
q0 = edge.getP0;
q1 = edge.getP1;
p0 = vertices(1).coord;
p = p0;
area2 = squaredArea(q0,q1,p0);
id = 1;
v = vertices(id);
volumeMin = 0;
for i = 2:length(vertices)
    pi = vertices(i).coord;
    volume = signedVolume(q0,q1,p,pi)
    if volume < 0
        p = pi;
        v = vertices(i);
        area2 = squaredArea(q0,q1,p);
    elseif (volume == 0)
        area2_ = squaredArea(q0,q1,pi);
        if area2_ > area2
            p = pi;
            v = vertices(i);
            area2 = area2_;
        end
    end
end
end

function area2 = squaredArea(p0,p1,p2)
u = p1 - p0;
v = p2 - p0;
area = cross(u,v)*0.5;
area2 = dot(area,area);
end
function volume = signedVolume(p0,p1,p2,p3)
u = p1 - p0;
v = p2 - p0;
z = p3 - p0;
normu = norm(u); normv = norm(v); normz = norm(z);
if ((normu == 0) || (normv == 0 ) || normz == 0)
    volume = 0.0;
    return ;
end
u = u / norm(u);
v = v / norm(v);
z = z / norm(z);
volume = det([u;v;z])
if abs(volume) < 10^-9
    volume = 0;
end
end
function m = showIds(triangles)
m = zeros(length(triangles),3);
for i = 1:length(triangles)
    t = triangles(i);
    m(i,:) = t.getIds();
end
m
end
function points = pointsInSphere()
n = 100;
teta = linspace(0,2*pi,n);
fi = linspace(0,pi,n);
x = cos(teta).*sin(fi);
y = sin(teta).*sin(fi);
z = cos(fi);
points = [x',y',z'];
theta=linspace(0,2*pi,40);
phi=linspace(0,pi,40);
[theta,phi]=meshgrid(theta,phi);
rho=1;
x=rho*sin(phi).*cos(theta);
y=rho*sin(phi).*sin(theta);
z=rho*cos(phi);
mesh(x,y,z)
end