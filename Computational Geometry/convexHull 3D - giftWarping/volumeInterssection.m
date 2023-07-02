clc;
clear;
close all;
A = 100*[0.2614    0.8959    0.4597;...
    0.2320    0.1959    0.8497;...
    0.7416    0.6622    0.1071;...
    0.8362    0.4549    0.6506;...
    0.9044    0.1107    0.6038;...
    0.0164    0.0912    0.9057;...
    0.5476    0.6442    0.2359;...
    0.6206    0.0790    0.6458];
B = 100*[ 0.1763    0.6189    0.9703;...
    0.6692    0.3318    0.4579;...
    0.0249    0.5231    0.9892;...
    0.7234    0.3969    0.5907;...
    0.5732    0.5077    0.1455;...
    0.2909    0.5011    0.1419;...
    0.4064    0.5931    0.8094;...
    0.4706    0.6741    0.7572];
filename = 'interseção entre cubos.gif';
%Default configurations
fig = figure;
hold on
angleInit = 30;
view(angleInit,30);
axis equal;
set(gca,'visible','off');
set(gcf,'color','white');
for i = 1:8
    plot3(A(i,1),A(i,2),A(i,3),'o','color','g');
end
pause(0.01);
count = 1;
frame = getframe(fig);
im{count} = frame2im(frame);
[fRot,count,im] = rot(angleInit,count,im,fig);
for i = 1:8
    plot3(B(i,1),B(i,2),B(i,3),'o','color','r');
end
pause(0.01);
[fRot,count,im] = rot(fRot,count,im,fig);

%setA = makeRand(8);
%setB = makeRand(8);
setA = Verts.empty;
setB = Verts.empty;
for i = 1:8
    setA(i) = Verts(i,A(i,:));
    setB(i) = Verts(i,B(i,:));
end
convexA = convexHull(setA);
convexB = convexHull(setB);
polyA = convexA.getConvexHull();
polyB = convexB.getConvexHull();
polyB.show([1,0,0]);
[fRot,count,im] = rot(fRot,count,im,fig);
[k1,av1] = convhull(A);
h2 = trisurf(k1,A(:,1),A(:,2),A(:,3),'FaceColor','green');
[fRot,count,im] = rot(fRot,count,im,fig);
%polyB.show([0,1,0]);
insideA = Verts.empty;
edges = polyB.getEdges;
i = 1;
insidePoints = Verts.empty;
nE = 8 + length(polyB.triangles) - 2;
coords = [];
for edge = edges
    [pts,i] = polyA.intersectionWithEdge(edge,i);
    n = length(pts);
    if n > 0
        insidePoints(end+1:end+n) = pts;
        for i = 1:n
            coords(end+1,1:3) = pts(i).coord;
        end
    end
end
for i = 1:length(insidePoints)
    plot3(coords(i,1),coords(i,2),coords(i,3),'o','color','b');
end
[fRot,count,im] = rot(fRot,count,im,fig);
delete(h2);
polyA.show([0,1,0]);
[fRot,count,im] = rot(fRot,count,im,fig);
[k1,av1] = convhull(coords);
trisurf(k1,coords(:,1),coords(:,2),coords(:,3),'FaceColor','cyan');
[fRot,count,im] = rot(fRot,count,im,fig);
convexC = convexHull(insidePoints);
polyC = convexC.getConvexHull();
% polyC.show('yellow');
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
    elseif (idx > 1) && (idx < length(im))        
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
    end
end

keyboard;
for j = 1:length(insideA)
    pt = insideA(j);
    showEdges(pt.edges,[0,0,1]);
    edgesFromPt = pt.edges;
    for edge = edgesFromPt
        edgeCoord = [pt.coord;edge.getIdOposto(pt).coord];
        pt = polyA.intersectionWithEdge(edgeCoord);
        if (length(pt) ~= 0)
            insideA(end+1) = pt;
        end
    end
end
convexC = convexHull(insideA);
polyC = convexC.getConvexHull();

function points = makeRand(n)
points(n) = Verts();
for i = 1:n
    points(i) = Verts(i,[rand,rand,rand]);
end
end
function bool = isInsideA(A,pt_)
bool = false;
for pt = A
    if pt == pt_
        bool = true;
        break;
    end
end
end
function showEdges(edges,color)
for edge = edges
    edge.show(color);
end
end
function [fRot,count,im] = rot(initangle,count,im,fig)
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
for i = 1:2:180
    fRot = initangle + i;
    view(fRot,30);
    pause(0.05);
    count = count + 1;
    frame = getframe(fig);
    im{count} = frame2im(frame);
end
end