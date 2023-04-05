classdef Edge < handle
    properties
        isSplited = false;
        isSplitedToRefine = false;
        hed1 = Hed.empty;
        hed2 = Hed.empty
        childrenId = zeros(1,2);
        id;
        mid;
        pRefined; %new points to refine
    end
    methods
        function this = Edge(hed1,hed2,id)
            if nargin == 3
                this.hed1 = hed1;
                this.hed2 = hed2;
                this.id = id;
            end
        end
        function out = getTwin(this,hedId)
            if this.hed1.id == hedId
                out = this.hed2;
            else
                out = this.hed1;
            end
        end
        function leafInsideNodes = getLeafInsideNodes(this,leafInsideNodes,edges)
            if (this.isSplited)
                leafInsideNodes(end+1) = this.mid;
                this.mid.markTemp = true;
                for i = 1:2
                    edgeChild = edges(this.childrenId(i));
                    leafInsideNodes = edgeChild.getLeafInsideNodes(leafInsideNodes,edges);
                end
            else
                return;
            end            
        end
        function plot(this,solid)
            p0  = solid.points(this.hed1.inc(1)).coord;
            p1  = solid.points(this.hed1.inc(2)).coord;
            line([p0(1),p1(1)], [p0(2),p1(2)], [p0(3),p1(3)], 'Color','green','linewidth',3);
        end
    end
end
