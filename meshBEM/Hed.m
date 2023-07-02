classdef Hed < handle
    properties
        isSplited = false;
        inc;
        id;
        children = Hed.empty;
        mid = Point.empty;
        
    end
    methods
        function this = Hed(inc,id)
            if nargin > 0
                this.inc = inc;
                this.id = id;
            end
        end
        function flip(this)
            this.inc = [ this.inc(2), this.inc(1) ] ;
        end
    end
end