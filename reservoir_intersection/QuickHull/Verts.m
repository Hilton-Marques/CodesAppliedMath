classdef Verts < handle
    properties
        id;
        coord;
        edges = Edge.empty;
        Marked = false;
    end
    methods
        function this = Verts(id,coord)
            if nargin > 0
                this.id = id;
                this.coord = coord;
            end
        end
        function show(this,color,size)
            if nargin == 2
                size = 5;
            end               
            plot3(this.coord(1),this.coord(2),this.coord(3),'o','MarkerFaceColor',color,'color','black','MarkerSize',size);
        end
    end
end