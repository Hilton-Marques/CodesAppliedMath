classdef Point < handle
    properties
        coord
        id
        hed0
        heds = Hed.empty;
    end
    methods
        function this = Point(coord,id)
            if nargin == 2
                this.coord = coord;
                this.id = id;
            end
        end
%         function hed = getHed(this,v)
%             hed = this.hed0;
%             v2 = hed.v2;
%             while (v2.id ~= v.id)
%                 hed = hed.getTwin().heNext
%                 v2 = hed.v2;
%             end
%         end
        function collectHedsNeigh(this)
            hed = this.hed0;
            firstId = hed.v2.id;
            this.heds(end+1) = hed;
            hed = hed.getTwin().heNext;
            v2 = hed.v2;
            while (v2.id ~= firstId)
                this.heds(end+1) = hed;
                hed = hed.getTwin().heNext;
                v2 = hed.v2;
            end
        end
        function hed = getHed(this,v)
            for i = 1: length(this.heds)
                hed = this.heds(i);
                if (hed.v2.id == v.id)
                    break;
                end
            end
        end
    end
end