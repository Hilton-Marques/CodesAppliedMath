classdef Face < handle
    properties
        hedInc;
        pointsInc;
        heds = Hed.empty;
        adjacentFaces;
        id;
        refineT6 = true;
        refinedPointsInc;
        q;
    end
    methods
        function this = Face(hedInc,heds,q)
            if nargin > 0
                this.hedInc = hedInc;
                this.getIndcPoints(heds);
                this.id = heds(hedInc).elId;
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
        function findAdjacents(this,heds,points,edges)
            elIds = [];
            id = this.id;
            for i = 1:length(this.heds)
                hed = this.heds(i);
                p2Target = hed.inc(2);
                hed = heds(edges(hed.edgeId).getTwin(hed.id).heNext);
                while true
                    p2 = hed.inc(2);
                    if p2 == p2Target
                        break;
                    end
                    newId = hed.elId;
                    if ~ismember(hed.elId,elIds)
                        elIds(end+1) = newId;
                    end
                    hed = heds(edges(hed.edgeId).getTwin(hed.id).heNext);
                end
            end
            this.adjacentFaces = elIds;
        end
        % refine for T6 and T9
        function newPoints = refine(this,edges,points)
            n = length(points);
            newPoints = Point.empty;
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
        function plot(this,points,color)
            if nargin == 2
                color = [0,0,1];
            end
            p = zeros(length(this.pointsInc),3);
            for i = 1:length(this.pointsInc)
                p(i,:) = points(this.pointsInc(i)).coord;
            end
            line([p(1,1),p(2,1)],[p(1,2),p(2,2)],...
                [p(1,3),p(2,3)],'Color', color);
            line([p(2,1),p(3,1)],[p(2,2),p(3,2)],...
                [p(2,3),p(3,3)],'Color',color);
            line([p(3,1),p(1,1)],[p(3,2),p(1,2)],...
                [p(3,3),p(1,3)],'Color',color);
            for i = 1:length(this.pointsInc)
                pi = p(i,:);
                plot3(pi(1),pi(2),pi(3),'o','Color','red');
            end
            if this.refineT6
            for i = 1:length(this.refinedPointsInc)
                pi = points(this.refinedPointsInc(i)).coord;
                plot3(pi(1),pi(2),pi(3),'o','Color','cyan');
            end
            end
        end
        function plotAdjacent(this,points,elements)
            n = length(this.adjacentFaces);
            for i = 1:n
                el = elements(this.adjacentFaces(i));
                el.plot(points,[1,0,0]);
            end
        end
    end
end