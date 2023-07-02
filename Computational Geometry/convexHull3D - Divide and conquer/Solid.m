classdef Solid < handle
    properties
        v
        tris = Face.empty;
        centroide = [0,0,0];
    end
    methods
        function this = Solid(pts)
            if nargin > 0
            this.v = pts;
            this.createSolid();
            this.initVerts();
            end
        end
        function createSolid(this)
            % triangle 1
            hed1 = Hed(this.v(1),this.v(2));
            edge1 = Edge(hed1);
            hed2 = Hed(this.v(2),this.v(3));
            edge2 = Edge(hed2);
            hed3 = Hed(this.v(3),this.v(1));
            edge3 = Edge(hed3);
            heds = [hed1,hed2,hed3];
            tri1 = Face(heds);
            this.v(1).hed0 = hed1;
            this.v(2).hed0 = hed2;
            this.v(3).hed0 = hed3;
            % triangle 2
            hed4 = edge1.getTwin(hed1);
            hed5 = Hed(hed4.v2, this.v(4));
            edge4 = Edge(hed5);
            hed6 = Hed(this.v(4), hed4.v1);
            edge5 = Edge(hed6);
            heds = [hed4,hed5,hed6];
            tri2 = Face(heds);
            this.v(4).hed0 = hed6;
            % triangle 3
            hed7 = edge2.getTwin(hed2);
            hed8 = edge5.getTwin(hed6);
            hed9 = Hed(this.v(4),hed7.v1);
            edge6 = Edge(hed9);
            heds = [hed7,hed8,hed9];
            tri3 = Face(heds);
            % triangle 4
            hed10 = edge3.getTwin(hed3);
            hed11 = edge6.getTwin(hed9);
            hed12 = edge4.getTwin(hed5);
            heds = [hed10,hed11,hed12];
            tri4 = Face(heds);                  
            this.tris = [tri1,tri2,tri3,tri4];
            
        end
        function initVerts(this)
            for pt = this.v
                pt.collectHedsNeigh();
                this.centroide = this.centroide + pt.coord;
            end
            this.centroide = this.centroide/length(this.v);
        end
        function show(this,color)
            if nargin == 1
                color = 'b';
            end
            for i = 1:length(this.tris)
                this.tris(i).show(color);
            end
        end
        function pt = getHigherPoint(this)
            zMax = this.v(1).coord(3);
            pt = this.v(1);
            for i = 2:4
                zTemp = this.v(i).coord(3);
                if (zTemp > zMax)
                    zMax = zTemp;
                    pt = this.v(i);
                end
            end
        end
    end
end