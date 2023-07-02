classdef convexHull < handle
    properties
        verts;
    end
    methods
        function this = convexHull(verts)
            this.verts = verts;            
        end
        function polyhedra = getConvexHull(this)
            vertices = this.verts;
            edges = Edge.empty;
            supportEdges = Edge.empty;
            triangles = Triangle.empty;
            triangle = this.findTriangleOnHull();
            supportEdges(end+1) = triangle.edges(1).getTwin();            
            for i = 2:3
                edgei = triangle.edges(i).getTwin();
                edges(end+1) = edgei;
                supportEdges(end+1) = edgei;
            end            
            triangles(end+1) = triangle;
            nel = 1;            
            while ( length(edges) ~= 0)
                edgei = edges(end);
                edges(end) = [];
                if (~edgei.isProcessed())
                    q = this.pivotOnEdge(edgei);
                    t = Triangle(edgei,q);                    
                    triangles(end+1) = t;
                    nel = nel + 1;
                    t.id = nel;                    
                    for i = 2:3
                        flag = false;                        
                        for k = 1:length(supportEdges)
                            edgek = supportEdges(k);
                            if (    (t.edges(i).p1.id == edgek.p0.id) && ...
                                    (t.edges(i).p0.id == edgek.p1.id) || ...
                                    (t.edges(i).p0.id == edgek.p0.id) && ...
                                    (t.edges(i).p1.id == edgek.p1.id)  )
                                edgek.processed = true;
                                flag = true;
                                t.edges(i) = edgek.twin;
                                break;
                            end
                        end
                        if ~flag
                            edgej = t.edges(i).getTwin;
                            edges(end+1) = edgej;
                            supportEdges(end+1) = edgej;
                        end
                    end
                    edgei.processed = true;
                end
            end
            polyhedra = Polyhedra(triangles);
        end
        function t = findTriangleOnHull(this)
            edge = this.findEdgeOnHull();
            p = this.pivotOnEdge(edge);
            t = Triangle(edge, p);
        end
        function v = pivotOnEdge(this,edge)
            vertices = this.verts;
            q0 = edge.getP0;
            q1 = edge.getP1;
            p0 = vertices(1).coord;
            p = p0;
            area2 = this.squaredArea(q0,q1,p0);
            id = 1;
            v = vertices(id);
            for i = 2:length(vertices)
                pi = vertices(i).coord;
                volume = this.signedVolume(q0,q1,p,pi);
                if volume < 0
                    p = pi;
                    v = vertices(i);
                elseif (volume == 0)
                    area2_ = this.squaredArea(q0,q1,pi);
                    if area2_ > area2
                        p = pi;
                        v = vertices(i);
                        area2 = area2_;
                    end
                end
            end
        end
        function edge = findEdgeOnHull(this)
            vertices = this.verts;
            pRight = vertices(1).coord;
            x = pRight(1);
            id = 1;
            for i =  2:length(vertices)
                pi = vertices(i).coord;
                x_ = pi(1);
                if (x_ > x)
                    x = x_;
                    pRight = pi;
                    id = vertices(i).id;
                end
            end
            suportVertex = Verts(length(vertices) + 1, pRight + [0,1,0]);
            edge = Edge ( suportVertex, vertices(id));
            r = this.pivotOnEdge(edge);
            edge = Edge(vertices(id),r);
        end
    end
    methods (Static)
        
        function area2 = squaredArea(p0,p1,p2)
            u = p1 - p0;
            v = p2 - p0;
            area = cross(u,v);
            area2 = dot(area,area);
        end
        function volume = signedVolume(p0,p1,p2,p3)
            u = p1 - p0;
            v = p2 - p0;
            z = p3 - p0;
            volume = det([u;v;z]);
            if abs(volume) < 0.0001
                volume = 0;
            end
        end
    end
end