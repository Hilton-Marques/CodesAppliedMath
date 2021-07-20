classdef Edge < handle
    properties
        isSplited = false;
        isSplitedToRefine = false;
        hed1 = Hed.empty;
        hed2 = Hed.empty
        childrenId = zeros(1,2);
        id;
        mid;
        pRefined; %new points to refine
    end
    methods
        function this = Edge(hed1,hed2,id)
            if nargin == 3
                this.hed1 = hed1;
                this.hed2 = hed2;
                this.id = id;
            end
        end
        function out = getTwin(this,hedId)
            if this.hed1.id == hedId
                out = this.hed2;
            else
                out = this.hed1;
            end
        end
    end
end
