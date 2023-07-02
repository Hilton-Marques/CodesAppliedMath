classdef Edge < handle
    properties
        hed1 = Hed.empty;
        hed2 = Hed.empty
    end
    methods
        function this = Edge(hed1)
            this.hed1 = hed1;
            this.hed2 = this.twin(hed1);
        end
        function out = twin(this)
            flipedInc = [this.hed1.inc(2),this.hed1.inc(1)];
            out = Hed(flipedInc,this.hed1.id);
            out.mid = this.mid;
        end
    end
end
