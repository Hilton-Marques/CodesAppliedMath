classdef Field < handle
    properties
        lmax;
        boundary;
        cells = cell.empty;
        levels = cell.empty;
        count = 1;
    end
    methods
        function this = Field(lmax,boundary)
            this.lmax = lmax;
            this.boundary = boundary;   
            q = boundary(1).nChild;
            nel = length(boundary)*(((q^lmax)-1)/(q-1));
            this.cells = cell(nel,2);
            this.levels = zeros(nel,lmax);
            this.start();
        end
        function start(this)
            for i = 1:length(this.boundary)
                bd = this.boundary(i);
                this.cells(this.count,1:2) = {bd,bd.key};
                this.levels(this.count,bd.level) = 1;
                this.count = this.count + 1;
                this.subdivide(bd);
            end
        end
        function subdivide(this,bd)
            if bd.level == this.lmax
                return
            end
            bd.subdivide();
            for i = 1:bd.nChild
                child = bd.child{i};
                this.cells{this.count} = child;
                this.levels(this.count,child.level) = 1;
                this.count = this.count + 1;
                this.subdivide(child);
            end
        end
    end
    
end