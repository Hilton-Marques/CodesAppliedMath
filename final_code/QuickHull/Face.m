classdef Face < handle
    properties
        hed1
        heds = HalfEdge.empty;
        pts = zeros(3,3);
        normal;
        eye_verts;
        pointsOnPositiveSide = Verts.empty;
        most_distant_point_dist = -1;
        most_distant_point = Verts.empty;
        far_point
        isVisibleFaceOnCurrentIter = false;
        was_visited_iter = -1;
        disabled = false;
        color = 'yellow'
        inc = zeros(1,3);
    end
    methods
        function this = Face(hed1)
            if nargin > 0
                this.hed1 = hed1;
                this.buildPts();
            end
        end
        function ids = getIds(this)
            ids = this.inc;
        end
        function buildPts(this)
            this.pts(1,:) = this.hed1.begin_vert.coord;
            hed = this.hed1;
            this.heds(end+1) = hed;
            this.inc(1) = hed.begin_vert.id;
            hed.face = this;
            for i = 2:3
                hed = hed.next;
                this.pts(i,:) = hed.begin_vert.coord;
                this.heds(end+1) = hed;
                hed.face = this;
                this.inc(i) = hed.begin_vert.id;
            end
            u = this.pts(2,:) - this.pts(1,:);
            v = this.pts(3,:) - this.pts(1,:);
            this.normal = cross(u,v);
            this.normal = this.normal / norm(this.normal);
        end
        function show(this,color) 
            if nargin == 1
                color = 'cyan';
            end
            this.color = color;
            trisurf([1,2,3],this.pts(:,1),this.pts(:,2),this.pts(:,3),'facecolor',color);            
            centroide = sum(this.pts) / 3;
            quiver3(centroide(1),centroide(2),centroide(3),0.5*this.normal(1),0.5*this.normal(2),0.5*this.normal(3),'color',color);
            for v = this.pointsOnPositiveSide
                v.show(color);
            end
            if ~isempty(this.most_distant_point)
                this.most_distant_point.show(color,10);
            end
        end
        function bool = addPoint(this,v)
            bool = false;
            d_signed = this.getDistanceFromPoint(v.coord);
            d = abs(d_signed);
            u = v.coord - this.pts(1,:);
            u_unit = u/norm(u);
            bias = abs(dot(u_unit,this.normal));
            
            if (d_signed > 0 && bias > 1.0e-10 )
                this.pointsOnPositiveSide(end+1) = v;
                if d > this.most_distant_point_dist
                    this.most_distant_point_dist = d;
                    this.most_distant_point = v;
                end
                bool = true;
            end
        end
        function d = getDistanceFromPoint(this,v)
            u = v - this.pts(1,:);
            d = dot(u,this.normal);
        end
        function showSet(this,color)
            for v = this.pointsOnPositiveSide
                v.show(color);
            end
        end
        function showNeigh(this)
            hed = this.hed1;
            for i = 1:3
                hed.twin.face.show('cyan');
                hed = hed.next;
            end
        end
    end
end