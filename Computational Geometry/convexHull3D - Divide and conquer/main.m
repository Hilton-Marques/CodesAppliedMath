clc;
clear;
close all;

[pts,coords] = makeRand(12);
coords = [-0.4840,-0.1904,0.2527;...
          -0.4607,-0.0198,0.3276;...
          -0.3182,-0.4551,0.2352;...
          -0.1886,0.3863,0.3120;...
          -0.0662,-0.2266,0.0507;...
          0.0124,-0.0688,-0.0186;...
          0.0134,0.3090,-0.5998;...
          0.0489,-0.0906,-0.1734];
[k1,av1] = convhull(coords(1:8,:));
[~,I] = sort(coords(:,1));
coords = coords(I,:);
pts = pts(I);
H1 = Solid(pts(1:4));
H2 = Solid(pts(5:8));
H3 = Solid(pts(9:12));
figure
hold on
view(30,30);
[k1,av1] = convhull(coords(1:8,:));
trisurf(k1,coords(1:8,1),coords(1:8,2),coords(1:8,3),'FaceColor','cyan');
axis equal;
hold off
figure
hold on
for i = 1:12
    plot3(coords(i,1),coords(i,2),coords(i,3),'o');
end
view(30,30);
H1.show();
H2.show();

% H3.show();

H4 = merge(H1,H2);
% for i = 1:4
%     tri1 = H1.tris(i);
%     tri2 = H2.tris(i);
%     if tri1.markTemp
%         tri1.show('black');
%     end
%     if tri2.markTemp
%         tri2.show('green');
%     end
% end
H4.show('cyan');


keyboard;

function [pts,coords] = makeRand(n)
pts = Point(n);
coords = zeros(n,3);
for i = 1:n
    coord = [rand,rand,rand];
    coords(i,:) = coord;
end
center = sum(coords,1)/n;
for i = 1:n
    coord = coords(i,:) - center;
    pts(i) = Point(coord,i);
    coords(i,:) = coord; 
end
end
function out = merge(H1,H2)
out = Solid();
%hedi = createFirstHed(H1,H2);
p1 = H1.getHigherPoint();
p2 = H2.getHigherPoint();
hedi = Hed(p1,p2);
Edge(hedi);
[newHed,tri] = pivotOnHed(H1,H2,hedi,true);
out.tris(end+1) = tri;
newHed = newHed.getTwin();
while ~(newHed.v1 == hedi.v1 && newHed.v2 == hedi.v2)
    [newHed,tri] = pivotOnHed(H1,H2,newHed,false);
    newHed = newHed.getTwin();
    out.tris(end+1) = tri;
end
for tri = H1.tris
    if ~(tri.markTemp)
        out.tris(end+1) = tri;
    end
end
for tri = H2.tris
    if ~(tri.markTemp)
        out.tris(end+1) = tri;
    end
end
end
function [newHed,tri] = pivotOnHed(H1,H2,hed,flag)
pts = [H1.v,H2.v];
q0 = hed.v1.coord;
q1 = hed.v2.coord;
p0 = pts(1).coord;
p = p0;
area2 = squaredArea(q0,q1,p0);
v = pts(1);
id = 1;
for i = 2:8
    pi = pts(i).coord;
    volume = signedVolume(q0,q1,p,pi);
    if volume < 0
        p = pi;
        v = pts(i);
        id = i;
    elseif (volume == 0)
        area2_ = squaredArea(q0,q1,pi);
        if area2_ > area2
            p = pi;
            v = pts(i);
            area2 = area2_;
            id = i;
        end
    end
end
if id > 4
    newHed = Hed(v,hed.v1);
    newEdge = Edge(newHed);
    oldHed = Hed(hed.v2,v);
    oldHed = hed.v2.getHed(v);
    if (flag)
        p1 = hed.v1.coord;
        p2 = hed.v2.coord;
        p3 = v.coord;
        u = p2 - p1;
        v = p3 - p1;
        normal = cross(u,v)/norm(cross(u,v));
        centroide = (p1+p2+p3)/3 ;
        %quiver3(centroide(1),centroide(2),centroide(3),normal(1),normal(2),normal(3));
        %quiver3(0,0,0,centroide(1),centroide(2),centroide(3));
%         if (dot(normal,centroide) < 0)
%             hed = hed.getTwin();
%             oldHed = oldHed.getTwin();
%             newHed = newHed.getTwin();
%         end
    end
    oldHed.face.markTemp = true;
    oldHed.face.show('black');
    %tri = Face([hed,oldHed,newHed]);
    tri = Face([hed.getTwin(),newHed.getTwin(),oldHed.getTwin()]);
    
else
    oldHed = Hed(v,hed.v1);
    oldHed = v.getHed(hed.v1);
    newHed = Hed(hed.v2,v);
    newEdge = Edge(newHed);
    if (flag)
        p1 = hed.v1.coord;
        p2 = hed.v2.coord;
        p3 = v.coord;
        u = p2 - p1;
        v = p3 - p1;
        normal = cross(u,v)/norm(cross(u,v));
        centroide = (p1+p2+p3)/3 ;
        centroideS = (H1.centroide + H2.centroide)*0.5;
        center = centroide - centroideS;
        %quiver3(centroide(1),centroide(2),centroide(3),normal(1),normal(2),normal(3));
        %quiver3(0,0,0,centroide(1),centroide(2),centroide(3));
%         if (dot(normal,center) < 0)
%             hed = hed.getTwin();
%             oldHed = oldHed.getTwin();
%             newHed = newHed.getTwin();
%         end
    end
    oldHed.show();
    oldHed.face.markTemp = true;
    oldHed.face.show('black');
    %tri = Face([hed,newHed,oldHed]);
    tri = Face([hed.getTwin(),oldHed.getTwin(),newHed.getTwin()]);
end
view(30,30);
tri.show('r');
tri.showNormal('g');
end

function v = pivotOnHedV(ptsH1,ptsH2,hed)
pts = [ptsH1,ptsH2];
q0 = hed.v1.coord;
q1 = hed.v2.coord;
p0 = pts(1).coord;
p = p0;
area2 = squaredArea(q0,q1,p0);
v = pts(1);
id = 1;
for i = 2:8
    pi = pts(i).coord;
    volume = signedVolume(q0,q1,p,pi);
    if volume < 0
        p = pi;
        v = pts(i);
        id = i;
    elseif (volume == 0)
        area2_ = squaredArea(q0,q1,pi);
        if area2_ > area2
            p = pi;
            v = pts(i);
            area2 = area2_;
            id = i;
        end
    end
end

end
function hed = createFirstHed(H1,H2)
p1 = H1.getHigherPoint();
p2 = H2.getHigherPoint();
hed = Hed(p1,p2);
p3 = pivotOnHedV(H1.v,H2.v,hed);
u = p2.coord - p1.coord;
v = p3.coord - p1.coord;
normal = cross(u,v);
centroide = (p1.coord+p2.coord+p3.coord)/3;
if (dot(normal,centroide) < 0 )
    hed = Hed(p2,p1);
end
Edge(hed);
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
if norm(u) ~= 0
    u = u/norm(u);
end
if norm(v)~= 0
    v = v/norm(v);
end
if norm(z)~=0
    z = z/norm(z);
end
volume = det([u;v;z]);
if abs(volume) < 0.00001
    volume = 0;
end
end