classdef Quad < handle
    properties
        indices;
    end
    methods
        function this = Quad(indices)
            if nargin > 1
            indices = indices;
            end
        end
        function plot(this,pts)
            for i = 1: 4
                pi = pts(this.indices(i),:);
                pj = pts(this.indices( mod(i,4)+1),:);
                line([pi(1),pj(1)],[pi(2),pj(2)], [pi(3),pj(3)], 'Color','red');
                
                
            end
        end
    end
end