classdef Face < handle
    properties
        hedInc;
        pointsInc;
        %points;
        heds = Hed.empty;
        adjacentFaces;
        id;
        typeEl;
        refinedPointsInc;
        %refinedPoints;
        children;
        father;
        q = []; %boundary condition
        bd = false;
        MEH;
        MEG;
        geometry;
        markTemp = false;
        colletPoints = [];
        adjacentElementsEdges = Face.empty;
        concave_adjacency = Face.empty;
        
    end
    methods
        function this = Face(hedInc,heds,q,typeEl)
            if nargin > 0
                this.hedInc = hedInc;
                this.getIndcPoints(heds);
                this.id = heds(hedInc).elId;                
                this.typeEl = typeEl;
                if ~(isempty(q))
                    this.q = q;
                    this.bd = true;
                end
            end
        end
        function getIndcPoints(this,heds)
            hedsEl(3) = Hed;
            pointsInc = zeros(1,3);
            hed = heds(this.hedInc);
            p1 = hed.inc(1);
            p2 = hed.inc(2);
            hedsEl(1) = hed;
            pointsInc(1) = p1;
            c = 1;
            while p2 ~= p1
                hed = heds(hed.heNext);
                c = c + 1;
                hedsEl(c) = hed;
                p2 = hed.inc(2);
                pointsInc(c) = hed.inc(1);
            end
            this.heds = hedsEl;
            this.pointsInc = pointsInc;
        end
        function findAdjacents(this,heds,edges,elements,points)
            elIds = [];
            id = this.id;
            for i = 1:length(this.heds)                
                hed = this.heds(i);
                if nargin == 4
                this.adjacentElementsEdges(end+1) = elements(edges(hed.edgeId).getTwin(hed.id).elId);
                end
                p2Target = hed.inc(2);
                hed = heds(edges(hed.edgeId).getTwin(hed.id).heNext);
                while true
                    p2 = hed.inc(2);
                    %points(p2).plot();
                    if p2 == p2Target
                        break;
                    end
                    newId = hed.elId;
                    elIds(end+1) = newId;
                    hed = heds(edges(hed.edgeId).getTwin(hed.id).heNext);
                end
            end
            for i = 1:length(this.concave_adjacency)
               elIds(end+1) = this.concave_adjacency(i).id;
            end
            this.adjacentFaces = unique(elIds);
        end
        function newPoints = refine(this,edges,points)
            n = length(points);
            newPoints = Point.empty;
            switch this.typeEl
                case 'Const'
                    n = n + 1;
                    coord = (points(this.pointsInc(1)).coord + ...
                        points(this.pointsInc(2)).coord + ...
                        points(this.pointsInc(3)).coord) /3;
                    pMid = Point(coord,n);
                    newPoints(end+1) = pMid;
                    this.refinedPointsInc(end+1) = n;
                case 'T3'
                    %do nothing
                case 'T6'
                    for i = 1:length(this.heds)
                        hed = this.heds(i);
                        edge = edges(hed.edgeId);
                        if ~edge.isSplitedToRefine
                            p1 = points(hed.inc(1));
                            p2 = points(hed.inc(2));
                            coord = (p1.coord + p2.coord)/2;
                            n = n+1;
                            pMid = Point(coord,n);
                            newPoints(end+1) = pMid;
                            this.refinedPointsInc(end+1) = n;
                            edge.pRefined = n;
                            edge.isSplitedToRefine = true;
                        else
                            this.refinedPointsInc(end+1) = edge.pRefined;
                        end
                    end
            end
        end
        function initialize(this,points)
            pointsTopology = points(this.pointsInc);
            refinedPoints = points(this.refinedPointsInc);
            
            % Pass boundary condition to points
%             for pt = pointsTopology
%                 if isempty(pt.q)
%                     pt.q = this.q;
%                 end
%             end
%             for pt = refinedPoints
%                 if isempty(pt.q)
%                     pt.q = this.q;
%                 end
%             end
            % Initialize geometry atributes
            totalPoints = [pointsTopology, refinedPoints];
            switch this.typeEl
                case 'Const'
                    this.geometry = Const(refinedPoints,pointsTopology);
                case 'T3'
                    this.geometry = T3(pointsTopology,pointsTopology);   
                case 'T6'
                    this.geometry = T6(totalPoints,totalPoints); 
            end
        end
        function leafElements = getLeafElements(this,leafElements)
            if isempty(this.children)
                leafElements(end+1) = this;
            end
            for child = this.children
                if isempty(child.children)
                    %é uma folha
                    leafElements(end+1) = child;
                else
                    leafElements = child.getLeafElements(leafElements);
                end
            end
        end
        function leafInsideNodes = findBoundary(this,heds,edges,elements,points)
            % just work if the adjacent elements were marked previously 
            
             boundaryEdges = Edge.empty;
             leafInsideNodes = Point.empty;
             
             adjIndices = this.adjacentFaces;
             adjIndices(end+1) = this.id;
             adjacentElements = elements(adjIndices);
             adjacentElements(end).markTemp = true;
             % find boundary first element
             flag = false;
             for adj = adjacentElements
                 adjHed = heds(adj.hedInc); 
                 adjHeds = [adjHed, heds(adjHed.heNext), heds(heds(adjHed.heNext).heNext)];
                 for hed = adjHeds
                     hedTwin = edges(hed.edgeId).getTwin(hed.id);
                     hedElement = elements(hedTwin.elId);
                     if (~hedElement.markTemp)
                         firstHed = hed;
                         flag = true;
                         break;
                     end
                 end
                 if flag
                     break;
                 end
             end
             v0 = firstHed.inc(1);
             boundaryEdges(end+1) = edges(firstHed.edgeId);
             pt = points(v0);
             pt.markTemp = true;
             leafInsideNodes(end+1) = pt;
             v1 = firstHed.inc(2);
             while (v0 ~= v1)
                 secondHed = heds(firstHed.heNext);
                 element = elements(secondHed.elId);
                 while (element.markTemp)
                     firstHed = secondHed;
                     edge = edges(firstHed.edgeId);
                     twinHed = edge.getTwin(firstHed.id);
                     secondHed = heds(twinHed.heNext);
                     element = elements(secondHed.elId);
                 end
                 v1 = firstHed.inc(1);
                 boundaryEdges(end+1) = edges(firstHed.edgeId);
                 pt = points(v1);
                 pt.markTemp = true;
                 leafInsideNodes(end+1) = pt;                 
             end
             for edge = boundaryEdges
                 leafInsideNodes = edge.getLeafInsideNodes(leafInsideNodes,edges);
             end
             % unMark the element 
             adjacentElements(end).markTemp = false;
        end
        function collectSources(this,sources)
            leafElements = Face.empty;
            leafElements = this.getLeafElements(leafElements);
            ids = [sources(:).id];
            if isempty(leafElements)
                leafElements = this;
            end
            for leaf = leafElements
                leaf.colletPoints(end+1 : end + length(ids)) = ids;
            end
        end
        function h = plot(this,color,alpha,flag, reverse)
            if nargin == 1
                color = 'b';
                flag = false;
                alpha = 0.5;
                reverse = +1;
            end
            if nargin == 2
                flag = false;
                alpha = 0.5;
                reverse = +1;
            end
            if nargin == 3
                flag = false;
                reverse = +1;
            end
            if nargin == 4
                reverse = +1;
            end
            p = zeros(length(this.pointsInc),3);
            for i = 1:length(this.pointsInc)
                p(i,:) = this.geometry.geometryNodes(i).coord;
            end
            if (flag)
                for i = 1:3
                    n = this.geometry.getNormal(0,0);
                    n = n/norm(n);
                    p(i,:) = p(i,:) + reverse*0.01*n';
                end
            end
            X = p(:,1);
            Y = p(:,2);
            Z = p(:,3);
%             line([p(1,1),p(2,1)],[p(1,2),p(2,2)],...
%                 [p(1,3),p(2,3)],'Color', [1,1,1,0.1]);
%             line([p(2,1),p(3,1)],[p(2,2),p(3,2)],...
%                 [p(2,3),p(3,3)],'Color',[1,1,1,0.1]);
%             line([p(3,1),p(1,1)],[p(3,2),p(1,2)],...
%                 [p(3,3),p(1,3)],'Color',[1,1,1,0.1]);
            h = trisurf([1,2,3],p(:,1),p(:,2),p(:,3),'FaceColor',color,'edgealpha',0.5,'facealpha',alpha);
            %h = fill3(X,Y,Z,color,'facealpha',alpha);
%             if (flag)
% %             for pi = this.geometry.fieldNodes
% %                 pCoords = pi.coord;
% %                 plot3(pCoords(1),pCoords(2),pCoords(3),'o','Color','red');
% %             end
%             end
            
        end
        function findConcave(this,range,fac)
            radius = this.geometry.getRadius(fac);
            for i = 1:3
                pt = this.geometry.geometryNodes(i).coord;
                for elementj = range
                    dist = elementj.geometry.getYc() - pt;
                    d = dot(dist,dist);
                    if d < radius
                        this.concave_adjacency(end+1) = elementj;
                    end
                end
            end
            this.concave_adjacency = unique(this.concave_adjacency);
        end
        function plotAdjacent(this,points,elements)
            n = length(this.adjacentFaces);
            for i = 1:n
                el = elements(this.adjacentFaces(i));
                el.plot(points,[1,0,0]);
            end
        end
        function plotChildren(this)
            yc = this.geometry.getYc();
            plot3(yc(1),yc(2),yc(3),'or','Color','red','MarkerFaceColor','r','MarkerSize',5);
            leafElements = Face.empty();
            leafElements = this.getLeafElements(leafElements);
            if isempty(leafElements)
                for i = 1:3
                    yCl = this.geometry.geometryNodes(i).coord;
                    plot(yCl,'o','Color','cyan');
                    u = yc - yCl;
                    quiver3(yCl(1),yCl(2),yCl(3),u(1),u(2),u(3),'Color','magenta');
                end
            else
            for (child = leafElements)
                child.plot('cyan',false);
                yCl = child.geometry.getYc();
                u = yc - yCl;
                quiver3(yCl(1),yCl(2),yCl(3),u(1),u(2),u(3),'Color','magenta');
                pause(0.3);
            end
            end
        end
        function plotHed(this,solid)
            heds = solid.heds;
            points = solid.points;
            for hed = this.heds
                hed.plot(heds,points,'r');
            end
        end
        
    end
end