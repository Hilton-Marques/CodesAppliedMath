classdef Hed < handle
    properties
        v1;
        v2;
        edge;
        heNext;
        face;
    end
    methods
        function this = Hed(v1,v2)
            this.v1 = v1;
            this.v2 = v2;
        end
        function show(this)
            line([this.v1.coord(1),this.v2.coord(1)],...
                [this.v1.coord(2),this.v2.coord(2)],...
                [this.v1.coord(3),this.v2.coord(3)], 'Color', 'g');
        end
        function hed = getTwin(this)
            hed = this.edge.getTwin(this);
        end
    end
end