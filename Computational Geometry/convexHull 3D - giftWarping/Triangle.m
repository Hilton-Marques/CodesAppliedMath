classdef Triangle < handle
    properties
        edges = Edge.empty;
        pts = zeros(3);
        verts = Verts.empty;
        id;
    end
    methods
        function this = Triangle(edge,p2)
            if nargin > 0
                this.pts(1,:) = edge.getP0;
                this.pts(2,:) = edge.getP1;
                this.pts(3,:) = p2.coord;
                
                this.verts(end+1) = edge.p0;
                this.verts(end+1) = edge.p1;
                this.verts(end+1) = p2;
                
                this.edges(end+1) = edge;
                this.edges(end+1) = Edge(this.verts(2), p2);
                this.edges(end+1) = Edge(p2, this.verts(1));

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

        
    end
end