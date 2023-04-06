classdef Mesh < handle
    properties
        m_verts;
        m_faces;
        m_edges;
        m_halfedges;
    end
    methods
        function this = Mesh()
        end
        function addVertex(this, vertex)
            this.m_verts(end+1) = vertex;
        end
        function addEdge(this, edge)
            this.m_edges(end +1) = edge;
            hed_1 = edge.getHed();
            hed_2 = edge.getHedTwin();
            this.addHalfEdge(this, hed_1);
            this.addHalfEdge(this, hed_2);
            this.addFace(hed.getFace());
        end
        function addHalfEdge(this, edge)
            this.m_halfedge(end+1) = edge;
        end
        function addFace(this, face)
            this.m_face(end+1) = face;
        end
    end
end