classdef Solver < Drawer
    properties
        m_vertices
        m_faces
        m_mesh
        m_geom
        m_parent_bb
    end
    methods
        function this = Solver(horizon_vertices, horizon_faces)
            this = this@Drawer();
            this.m_vertices = horizon_vertices;
            this.m_faces = horizon_faces;
            this.m_mesh = ManifoldSurfaceMesh(horizon_faces);
            this.m_geom = Geometry(horizon_vertices);
            this.showMesh();
            this.m_parent_bb = this.getBB(this.m_vertices);
            this.setBB(this.m_parent_bb);
            this.showBoundaryLoops();
        end
        function showMesh(this)
            trisurf(this.m_faces, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
                'FaceAlpha',0.7,'EdgeColor','none','FaceColor', this.m_blue);
        end

        function focusCamera(this, vertices_ids)
            coords = this.m_geom.inputVertexPosition(vertices_ids);
            bb = this.getBB(coords);
            this.walkTrough(bb);
            this.walkTrough(this.m_parent_bb);
        end
    end
end