classdef Face < handle
    properties
        heds;
        markTemp = false;
    end
    methods
        function this = Face(heds)
            this.heds = heds;
            for i = 1:3
                this.heds(i).face = this;
                this.heds(i).heNext = heds(mod(i,3)+1);
            end
        end
        function show(this,color)
            if nargin == 1
                color = 'blue';
            end
            for i = 1:3
               line([this.heds(i).v1.coord(1),this.heds(i).v2.coord(1)],...
               [this.heds(i).v1.coord(2),this.heds(i).v2.coord(2)],...
               [this.heds(i).v1.coord(3),this.heds(i).v2.coord(3)],'Color',color);               
            end
            normal = this.normal();
            centroide = this.getCentroide();
            %quiver3(centroide(1),centroide(2),centroide(3),normal(1),normal(2),normal(3));
        end
        function normal = normal(this)
            u = this.heds(1).v2.coord - this.heds(1).v1.coord;
            v = this.heds(2).v2.coord - this.heds(2).v1.coord;
            normal = cross(u,v);
            normal = normal/norm(normal);
        end
        function  centroide = getCentroide(this)
            p1 = this.heds(1).v1.coord;
            p2 = this.heds(1).v2.coord;
            p3 = this.heds(2).v2.coord;
            centroide = (p1+p2+p3)/3;
        end
        function showNormal(this,color)
            centroide = this.getCentroide();
            normal = this.normal();         
            quiver3(centroide(1),centroide(2),centroide(3),normal(1),normal(2),normal(3),'Color',color);
        end
        function bool = correctOrientation(this)
            normal = this.normal();
            centroide = this.getCentroide();
            bool = false;
            if (dot(centroide,normal) > 0)
                bool = true;
            end
        end
    end
end