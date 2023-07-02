classdef HalfPlane < handle
    properties
        pts(4,1) = Vertex();
        normal;
        d;
    end
    methods
        function this = HalfPlane(pts, n, d)
            if nargin == 1
                this.pts = pts;
                this.normal = this.getNormal();
            end
            if nargin == 3
                this.normal = n;
                this.d = d;
            end
        end
        function n = getNormal(this)
            p1 = this.pts(2).coord - this.pts(1).coord;
            p2 = this.pts(3).coord - this.pts(1).coord;
            n = cross(p2,p1);
            n = n/norm(n);
        end
        function angle = getAngle(this)
            n = this.getNormal();
            i = [1,0,0];
            angle = acos(dot(n,i));
            if det([i;n;cross(n,i)]) < 0
                angle = 2*pi - angle;
            end
        end
        function translate(this,trans)
            for i = 1:4
                pt = this.pts(i);
                pt.coord = pt.coord - trans*this.normal;
            end
            
        end
        function out = getD(this)
            out = dot(this.normal,this.pts(1).coord);
        end
        function show(this,color)
            if nargin == 1
                color = 'blue';
            end
            for i = 1:4
                line([this.pts(i).coord(1),this.pts(mod(i,4)+1).coord(1)], ...
                    [this.pts(i).coord(2),this.pts(mod(i,4)+1).coord(2)], ...
                    [this.pts(i).coord(3),this.pts(mod(i,4)+1).coord(3)], 'Color',color);
            end
        end
        function pt = transform(this)
            d = dot(this.normal,this.pts(1).coord);
            pt = this.normal*(1/d);
        end
    end
end