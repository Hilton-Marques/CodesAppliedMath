classdef Edge < handle
    properties
        p0;
        p1;
        processed = false;
        twin;
    end
    methods
        function this = Edge(p0, p1)
            if nargin > 0
                this.p0 = p0;
                this.p1 = p1;
                this.startTwin();              
            end
        end
        function startTwin(this)
            edge = Edge();
            edge.p0 = this.p1;
            edge.p1 = this.p0;
            this.twin = edge;
            edge.twin = this;
        end
        function p0 = getP0(this)
            p0 = this.p0.coord;
        end
        function p1 = getP1(this)
            p1 = this.p1.coord;
        end
        function bool = isProcessed(this)
            bool = this.processed;
        end
        function edge = getTwin(this)
            edge = this.twin;
            %update points
            this.p0.edges(end+1) = edge;
            this.p1.edges(end+1) = edge;
        end
        function out = getIdOposto(this, pt)
            out = this.p1;
            if this.p1.id == pt.id
                out = this.p0;
            end
            
        end
        function show(this,color)
            if nargin == 1
                color = 'red';
            end
            line([this.p0.coord(1),this.p1.coord(1)],...
                [this.p0.coord(2),this.p1.coord(2)],...
                [this.p0.coord(3),this.p1.coord(3)],...
                'color',color);
        end
    end
end