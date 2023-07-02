classdef HalfEdge < handle
    properties
        begin_vert
        end_vert
        face_id
        next
        twin
        face
    end
    methods
        function this = HalfEdge()
            if nargin > 0

            end
        end
        function show(this,color)
            pi = this.begin_vert.coord;
            pj = this.end_vert.coord;
            line([pi(1),pj(1)],[pi(2),pj(2)],[pi(3),pj(3)],'color',color);
        end
    end
end