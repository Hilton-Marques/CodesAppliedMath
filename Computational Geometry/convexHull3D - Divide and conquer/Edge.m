classdef Edge < handle
    properties
        hed1;
        hed2;
    end
    methods
        function this = Edge(hed1)
            this.hed1 = hed1;
            this.hed2 = Hed(hed1.v2, hed1.v1);      
            this.hed1.edge = this;
            this.hed2.edge = this;
        end
        function out = getTwin(this,hed)
            if this.hed1 == hed
                out = this.hed2;
            else
                out = this.hed1;
            end
        end
    end
end
