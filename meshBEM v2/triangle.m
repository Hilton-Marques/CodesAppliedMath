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
heds(1) = Hed([1,2],1,1,1,2);
heds(2) = Hed([2,3],2,2,1,3);
heds(3) = Hed([3,1],3,3,1,1);
heds(4) = Hed([4,1],4,4,2,5);
heds(5) = Hed([1,3],5,3,2,6);
heds(6) = Hed([3,4],6,6,2,4);
heds(7) = Hed([3,2],7,2,3,8);
heds(8) = Hed([2,4],8,5,3,9);
heds(9) = Hed([4,3],9,6,3,7);
heds(10) = Hed([2,1],10,1,4,11);
heds(11) = Hed([1,4],11,4,4,12);
heds(12) = Hed([4,2],12,5,4,10);
elements(4) = Face;
elements(1) = Face(1,heds);
elements(2) = Face(4,heds);
elements(3) = Face(7,heds);
elements(4) = Face(10,heds);
edges(6) = Edge;
edges(1) = Edge(heds(1),heds(10),1);
edges(2) = Edge(heds(2),heds(9),2);
edges(3) = Edge(heds(3),heds(5),3);
edges(4) = Edge(heds(4),heds(11),4);
edges(5) = Edge(heds(8),heds(12),5);
edges(6) = Edge(heds(6),heds(9),6);
nL = 4; % 4 levels
allFaces = cell(nL,1);
allFaces{1} = elements;
nHeds = 12;
nPts = 4;
nEdges = 6;
nEl = 4;
shouldRefine = false;
for i = 2:nL
    Elements_m = allFaces{i-1};
    nElNv = length(Elements_m);
    newElements(4*nElNv) = Face;
    count = 1;
    for j = 1:nElNv
        el = Elements_m(j);
        localHeds = Hed.empty;
        localPoints = Point.empty;
        for k = 1:length(el.heds)
            hed_m = el.heds(k);
            edge_m = edges(hed_m.edgeId);
            if ~edge_m.isSplited
                coord = (points(hed_m.inc(1)).coord + points(hed_m.inc(2)).coord )/2;
                nPts = nPts + 1;
                pt = Point(coord,nPts);
                points(end+1) = pt;
                localPoints(end+1) = pt;
                for l = 1:2 % number of subdivisions
                    nEdges = nEdges + 1;
                    nHeds = nHeds + 1;
                    inc1 = [hed_m.inc(l),nPts];
                    inc2 = [nPts,hed_m.inc(l)];
                    incs = [inc1;inc2];
                    newHed = Hed(incs(l,:),nHeds,nEdges);
                    nHeds = nHeds + 1;
                    newHedTwin = Hed(incs(mod(l,2)+1,:),nHeds,nEdges);
                    newEdge = Edge(newHed,newHedTwin,nEdges);
                    localHeds(end+1) = newHed;
                    heds(end+1) = newHed;
                    heds(end+1) = newHedTwin;
                    edges(end+1) = newEdge;
                    edge_m.childrenId(l) = nEdges;
                end
                edge_m.isSplited = true;
                edge_m.mid = pt;
            else
                localPoints(end+1) = edge_m.mid;
                for l = 2:-1:1 % the heds reverse the order
                    edge_new = edges(edge_m.childrenId(l));
                    localHeds(end+1) = edge_new.hed2;
                end
            end
        end
        internalHeds(3) = Hed;
        for k = 1:3
            nEdges = nEdges + 1;
            nHeds = nHeds + 1;
            p2 = localPoints(mod(k,3)+1).id;
            p1 = localPoints(k).id;
            newHed = Hed([p1,p2],nHeds,nEdges);
            nHeds = nHeds + 1;
            newHedTwin = Hed([p2,p1],nHeds,nEdges);
            internalHeds(k) = newHed;
            heds(end+1) = newHed;
            heds(end+1) = newHedTwin;
            newEdge = Edge(newHed,newHedTwin,nEdges);
            edges(end+1) = newEdge;
        end
        %reorder
        internalHeds = internalHeds([3,1,2]) ;
        % Update hNext for border triangles
        c1 = 1;
        c2 = 5;
        for l = 1:2:5
            nEl = nEl + 1;
            firstId = localHeds(l).id;
            idMidle = edges(internalHeds(c1).edgeId).hed2.id;
            localHeds(l).heNext = idMidle;
            localHeds(l).elId = nEl;
            idBefore = mod(6 + l - 2,6)+1;
            heds(idMidle).heNext = localHeds(idBefore).id;
            heds(idMidle).elId = nEl;
            localHeds(mod(c2,6)+1).heNext = firstId;
            localHeds(mod(c2,6)+1).elId = nEl;
            c1 = c1 + 1;
            c2 = c2 + 2;
        end
        % Update hNext for internal triangle
        nEl = nEl + 1;
        for l = 1:3
            internalHeds(l).heNext = internalHeds(mod(l,3)+1).id;
            internalHeds(l).elId = nEl;
        end
        %Create new Elements
        elementsChildren(4) = Face;
        elementsChildren(1) = Face(localHeds(1).id,heds);
        elementsChildren(2) = Face(localHeds(3).id,heds);
        elementsChildren(3) = Face(localHeds(5).id,heds);
        elementsChildren(4) = Face(internalHeds(1).id,heds);
        newElements(count:count+3) = elementsChildren;
        count = count + 4;
    end
    if shouldRefine
        for k = 1:length(newElements)
            el = newElements(k);
            newPoints = el.refine(edges,points);
            % As new points was created there is a need to updtate nPts
            points(end+1:end+length(newPoints)) = newPoints;
            nPts = length(points);
        end
    end
    allFaces{i} = newElements;
end
% vectorize elements
elements(nEl) = Face;
count = 1;
for i = 1:length(allFaces)
    els = allFaces{i};
    elements(count:count+length(els)-1) = els;
    count = count+length(els);
end
%% Plot
nv = 4;
axis equal;
ax = gca;
set(gca,'XColor', 'none','YColor','none','ZColor','none')
view(30,70);
elementsNv  = allFaces{nv};
hold on
for i = 1:length(elementsNv)
    el = elementsNv(i);
    el.plot(points);
    if nv > 1
        el.findAdjacents(heds,points,edges);
    end
end
if nv > 1
    el = elementsNv(1);
    el.plotAdjacent(points,elements);
    el.plot(points,[0,1,0]);
end
keyboard

