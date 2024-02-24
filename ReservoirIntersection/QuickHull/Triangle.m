classdef Triangle < handle
    properties
        edges = Edge.empty;
        pts = zeros(3);
        verts = Verts.empty;
        id;
        v1;
        v2;
        v3;
        normal;
    end
    methods
        function this = Triangle(v1,v2,v3)
            if nargin > 0
                this.v1 = v1;
                this.v2 = v2;
                this.v3 = v3;
                this.pts = [this.v1.coord;this.v2.coord;this.v3.coord];
                u = this.v2.coord - this.v1.coord ;
                v = this.v3.coord - this.v1.coord;
                normal = cross(u,v);
                this.normal = normal / norm(normal);                
            end    
        end
        
        function show(this, color)
            if nargin == 1
                color = 'blue';
            end
            trisurf([1,2,3],this.pts(:,1),this.pts(:,2),this.pts(:,3),'FaceColor','cyan','Facealpha',0.8);
            for i = 1:3
                
                line([this.pts(i,1),this.pts(mod(i+1,3)+1,1)], ...
                     [this.pts(i,2),this.pts(mod(i+1,3)+1,2)], ...
                     [this.pts(i,3),this.pts(mod(i+1,3)+1,3)], 'Color',color);
            end
        end
        function ids = getIds(this)
            ids = [this.verts(1).id,this.verts(2).id,this.verts(3).id];            
        end
        function normal = getNormal(this)
            u = this.pts(2,:) - this.pts(1,:);
            v = this.pts(3,:) - this.pts(1,:);
            u = u/norm(u);
            v = v/norm(v);
            normal = cross(u,v);
        end
        function centroide = getCentroide(this)
            centroide = (this.pts(1,:) + this.pts(2,:) + this.pts(3,:))/3;
        end
        function d = getDistanceFromPoint(this,v)
            u = v - this.v1.coord;
            d = (dot(u,this.normal));
        end

        
    end
end