classdef Triangle < handle
    properties
        points = zeros(3,2);
        neigh = Triangle.empty;
        edges
        % Circuncircle properties
        center
        radius
    end
    methods
        function this = Triangle(p1,p2,p3,neigh)
            if nargin > 3
                this.neigh = neigh;
            end
            if nargin > 0
                this.points = [p1;p2;p3];
                this.createCircumcircle();
                edges(3) = Edge;
                edges(1) = Edge(p2,p3);
                edges(2) = Edge(p3,p1);
                edges(3) = Edge(p1,p2);
                this.edges = edges;
            end
        end
        function createCircumcircle(this)
            p1 = this.points(1,:);
            p2 = this.points(2,:);
            p3 = this.points(3,:);
            l1 = p2 - p1;
            l2 = p3 - p1;
            a = l2;
            b = l1;
            ac = [-a(2),a(1)];
            bc = [-b(2),b(1)];
            rc = (dot(b,b)*ac - dot(a,a)*bc)/(2*dot(ac,b));
            this.center = p1 + rc;
            this.radius = dot(rc,rc);
        end
        function out = isInCircumcircle(this,p)
            d = dot(p-this.center,p-this.center);
            out = false;
            if (d <= this.radius)
                out = true;
            end
        end
        function plot(this,c)
            for i = 1:3
                edge = this.edges(i);
                p1 = edge.p1;
                p2 = edge.p2;
                line([p1(1),p2(1)],[p1(2),p2(2)],'Color',c);
                pMed = (p1+p2)/2;
                %text(pMed(1),pMed(2),num2str(i));
            end
        end
    end
end