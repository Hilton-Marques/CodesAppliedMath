% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;
%% Domain
%Piramid
points(4) = Point;
points(1) = Point([-3;-1.5;0],1);
points(2) = Point([3;-1.5;0],2);
points(3) = Point([0;1.5;0],3);
points(4) = Point([0;0;4],4);

%% Split the domain
heds(6) = Hed; %number of edges
heds(1) = Hed([1,2],1);
heds(2) = Hed([2,3],2);
heds(3) = Hed([3,1],3);
heds(4) = Hed([4,1],4);
heds(5) = Hed([2,4],5);
heds(6) = Hed([3,4],6);
Elements(4) = Face;
Elements(1) = Face([1,2,3],heds);
Elements(2) = Face([4,3,6],heds);
Elements(3) = Face([2,5,6],heds);
Elements(4) = Face([1,4,5],heds);
nL = 4; % 4 levels
allFaces = cell(nL,1);
allFaces{1} = Elements;
nHed = 6;
nPts = 4;
for i = 2:nL
    Elements_m = allFaces{i-1};
    nEl = length(Elements_m);
    newElements(4*nEl) = Face;
    count = 1;
    for j = 1:nEl
        el = Elements_m(j);
        localHeds = Hed.empty;
        localPoints = Point.empty;
        bool = [];
        for k = 1:3
            hed_m = heds(el.hedInc(k));
            if ~hed_m.isSplited
                coord = (points(hed_m.inc(1)).coord + points(hed_m.inc(2)).coord )/2;
                nPts = nPts + 1;
                pt = Point(coord,nPts);
                nHed = nHed + 1;
                newHed1 = Hed([hed_m.inc(1),nPts],nHed);
                localHeds(end+1) = newHed1;
                nHed = nHed + 1;
                newHed2 = Hed([nPts,hed_m.inc(2)],nHed);
                localHeds(end+1) = newHed2;
                localPoints(end + 1) = pt;
                children(2) = Hed;
                children(1) = newHed1;
                children(2) = newHed2;
                hed_m.children = children;
                heds(end + 1:end+2) = children;
                points(end + 1) = pt;
                hed_m.mid = pt;
                hed_m.isSplited = true;
            else
                hed_m.children(2).flip;
                hed_m.children(1).flip
                localHeds(end+1) = hed_m.children(2);
                localHeds(end+1) = hed_m.children(1);
                localPoints(end+1) = hed_m.mid;
                bool(end+1) = hed_m.children(2).id;
                bool(end+1) = hed_m.children(1).id;
            end
        end
        for k = 1:3
            nHed = nHed + 1;
            p2 = localPoints(k).id;
            p1 = localPoints(mod(k,3)+1).id;
            newHed = Hed([p1,p2],nHed);
            localHeds(end+1) = newHed;
            heds(end+1) = newHed;
        end
        elementsChildren(4) = Face;
        elementsChildren(1) = Face([localHeds(1).id,localHeds(9).id,localHeds(6).id],heds);
        elementsChildren(2) = Face([localHeds(2).id,localHeds(3).id,localHeds(7).id],heds);
        elementsChildren(3) = Face([localHeds(4).id,localHeds(5).id,localHeds(8).id],heds);
        elementsChildren(4) = Face([localHeds(9).id,localHeds(7).id,localHeds(8).id],heds);
        newElements(count:count+3) = elementsChildren;
        count = count + 4;
        if ~isempty(bool)
            for k = 1:length(bool)
                heds(bool(k)).flip;
            end
        end
    end
    allFaces{i} = newElements;
end
%% Plot
nv = 3;
axis equal;
view(30,70);
elements  = allFaces{nv};
hold on
for i = 1:length(elements)
    el = elements(i);
    el.plot(points);
end


function plotTri(conec,points)
p = points(conec);
line([p(1).coord(1),p(2).coord(1)],[p(1).coord(2),p(2).coord(2)],...
    [p(1).coord(3),p(2).coord(3)]);
line([p(2).coord(1),p(3).coord(1)],[p(2).coord(2),p(3).coord(2)],...
    [p(2).coord(3),p(3).coord(3)]);
line([p(3).coord(1),p(1).coord(1)],[p(3).coord(2),p(1).coord(2)],...
    [p(3).coord(3),p(1).coord(3)]);
end
function out = flipHed(hed)
flipedInc = [hed.inc(2),hed.inc(1)];
out = Hed(flipedInc,hed.id);
out.mid = hed.mid;
end