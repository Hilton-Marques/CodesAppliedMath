classdef Mesh < handle
    properties
        points
    end
    methods
        function this = Mesh(points)
            this.points = points; 
            for i=1:size(points,1)
            end
        end
    end
end