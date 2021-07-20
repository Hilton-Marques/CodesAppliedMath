classdef Hed < handle
    properties
        inc;
        id;
        elId;
        edgeId;
        heNext;
    end
    methods
        function this = Hed(inc,id,edgeId,elId,heNext)
            if nargin == 5
                this.inc = inc;
                this.id = id;
                this.edgeId = edgeId;
                this.heNext = heNext;
                this.elId = elId;
            elseif nargin == 3
                this.inc = inc;
                this.id = id;
                this.edgeId = edgeId;
            end
        end
        function flip(this)
            this.inc = [ this.inc(2), this.inc(1) ] ;
        end
    end
end