classdef Vertex < handle
    properties 
        coord;
        id;
    end
    methods
        function this = Vertex(coord,id)
            if nargin > 0
                this.coord = coord;
                this.id = id;
            end
        end
    end
end