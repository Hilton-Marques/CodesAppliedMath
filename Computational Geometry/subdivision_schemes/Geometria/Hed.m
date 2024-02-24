classdef Hed < handle
    properties
        inc;
        id;
        elId;
        edgeId = -1;
        heNext;
    end
    methods
        function this = Hed(inc,id,elId,heNext,edgeId)
            if nargin == 5
                this.inc = inc;
                this.id = id;
                this.edgeId = edgeId;
                this.heNext = heNext;
                this.elId = elId;
            elseif nargin == 3
                this.inc = inc;
                this.id = id;
                this.edgeId = elId;
            elseif nargin == 4
                this.inc = inc;
                this.id = id;
                this.elId = elId;
                this.heNext = heNext;
            end
        end
        function flip(this)
            this.inc = [ this.inc(2), this.inc(1) ] ;
        end
        function plot(this,heds, points, color)
            
            p0  = points(this.inc(1)).coord;
            p1  = points(this.inc(2)).coord;
            u = p1 - p0;    
            L = norm(u);
            uc = u / L;
            heNext = heds(this.heNext);
            p2 = points(heNext.inc(2)).coord;
            v = p2 - p1;
            n = cross(u,v);
            t = cross(n,u);
            t = t / norm(t);
            p0 = 0.1*t + p0 + 0.2*u ;
            ud = 0.8 * L * uc;
            quiver3(p0(1),p0(2),p0(3),ud(1),ud(2),ud(3),'Color',color);
        end
    end
end