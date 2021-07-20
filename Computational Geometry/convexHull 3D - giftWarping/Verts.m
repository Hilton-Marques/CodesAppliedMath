classdef Verts < handle
    properties
        id;
        coord;
        edges = Edge.empty;
    end
    methods
        function this = Verts(id,coord)
            if nargin > 0
                this.id = id;
                this.coord = coord;
            end
        end
    end
end