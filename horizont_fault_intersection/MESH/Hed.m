classdef Hed < handle
    properties
        m_inc;
        m_id;
        m_face
        m_edge
        m_heNext;
        m_p1;
        m_p2;
    end
    methods
        function this = Hed(inc,id,p1,p2)
            this.m_inc = inc;
            this.m_id = id;
            this.m_p1 = p1;
            this.m_p2 = p2;
        end
        function setEdge(this, edge)
           this.m_edge = edge; 
        end
        function setFace(this, face)
            this.m_face = face;
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