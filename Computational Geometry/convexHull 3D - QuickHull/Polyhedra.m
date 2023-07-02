classdef Polyhedra < handle
    properties
        triangles
    end
    methods
        function this = Polyhedra(triangles)
            this.triangles = triangles;
        end
        function bool = isInside(this,pt)
            bool = true;
            for t = this.triangles
                volume = this.signedVolume(t.pts(1,:),t.pts(2,:), t.pts(3,:), pt.coord );
                if volume < 0
                    bool = false;
                    return;
                end
            end
        end
        function edges = getEdges(this)
            edges = Edge.empty;
            for t = this.triangles
                for i = 1:3
                    edges(end+1) = t.edges(i);
                end                    
            end
            edges = unique(edges);
        end
        function [pts,i] = intersectionWithEdge(this,edge,i)
            pts = Verts.empty;
            t = edge.getP0();
            v = edge.getP1() - edge.getP0();
            a = this.pluckCoordinates(edge.getP0(),edge.getP1());
            for tri = this.triangles
                l1 = this.pluckCoordinates(tri.pts(1,:),tri.pts(2,:));
                l2 = this.pluckCoordinates(tri.pts(2,:),tri.pts(3,:));
                l3 = this.pluckCoordinates(tri.pts(3,:),tri.pts(1,:));
                s1 = this.side(a,l1);
                s2 = this.side(a,l2);
                s3 = this.side(a,l3);
                sign1 = sign(s1);
                if sign(s2) ~= sign1 || sign(s3) ~= sign1
                    continue;
                end
                l2 = this.pluckCoordinates(edge.getP1(),tri.pts(1,:));
                l3 = this.pluckCoordinates(tri.pts(1,:),edge.getP0());
                l4 = this.pluckCoordinates(tri.pts(2,:),tri.pts(3,:));
                s1 = this.side(l4,l3);
                s2 = this.side(l4,l2);
                if sign(s1) ~= sign(s2)
                    continue;
                end
                n = tri.getNormal();
                d = dot(tri.pts(1,:),n);
                lam = (d - dot(n,t))/dot(n,v);
                if lam >= 0 && lam <= 1
                    coord = t + lam*v;
                    plot3(coord(1),coord(2),coord(3),'o');
                    newPt = Verts(i,coord);
                    pts(end+1) = newPt;
                    i = i + 1;
                end
            end
            
        end
        function show(this, color)
            for t = this.triangles
                t.show(color);
            end
        end
    end
    
    methods (Static)
        function volume = signedVolume(p0,p1,p2,p3)
            u = p1 - p0;
            v = p2 - p0;
            z = p3 - p0;
            volume = det([u;v;z]);
            if abs(volume) < 0.0001
                volume = 0;
            end
        end
        function out = pluckCoordinates(p,q)
            out = zeros(6,1);
            out(1) = p(1)*q(2) - q(1)*p(2);
            out(2) = p(1)*q(3) - q(1)*p(3);
            out(3) = p(1) - q(1);
            out(4) = p(2)*q(3) - q(2)*p(3);
            out(5) = p(3) - q(3);
            out(6) = q(2) - p(2);
        end
        function out = side(a,b)
            out = a(1)*b(5) + a(2)*b(6) + a(3)*b(4) + ...
                a(4)*b(3) + a(5)*b(1) + a(6)*b(2);
            if (abs(out) < 0.00001)
                out = 0;
            end
                
        end
    end
end