classdef Point < handle
    properties
        coord
        id
    end
    methods
        function this = Point(coord,id)
            if nargin == 2
                this.coord = coord;
                this.id = id;
            end
        end
    end
end