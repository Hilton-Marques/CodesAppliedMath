classdef Fault < handle
    properties
        m_vertices
        m_faces
        m_mesh
        m_geom
    end
    methods
        function this = Fault(fault_vertices, fault_faces)
            this.m_vertices = fault_vertices;
            this.m_faces = fault_faces;
            this.m_mesh = ManifoldSurfaceMesh(horizon_faces);
            this.m_geom = Geometry(fault_vertices);
        end
        function showMesh(this, color)
            trisurf(this.m_faces, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
                'FaceAlpha',0.7,'EdgeColor','none','FaceColor', color);
        end

        function showBoundaryLoops(this)
            n_loops = this.m_mesh.m_nboundary_loops;
            colors = rand(n_loops,3);

            if (n_loops > 0)
                %the boundary loops are at the end of faces array
                first_loop = this.m_mesh.m_nfaces + 1;

                for id_loop = first_loop:(first_loop + n_loops - 1)
                    first_hed = this.m_mesh.m_fHalfEdgeArr(id_loop);
                    prev_hed = first_hed;
                    next_hed = ManifoldSurfaceMesh.INVALID_IND;
                    %save vertices for bounding box
                    vertices_id = [];
                    color = colors(id_loop - this.m_mesh.m_nfaces, :);
                    while (next_hed ~= first_hed)
                        next_hed = this.m_mesh.m_heNextArr(prev_hed);
                        vertices_id = [vertices_id,...
                            this.m_mesh.m_heVertexArr(prev_hed)];
                        % plot edge and vertices
                        this.plotEdge(prev_hed,next_hed,color);
                        prev_hed = next_hed;
                    end
                    %                     for i = 1:10
                    %                         this.update();
                    %                     end
                    %this.focusCamera(vertices_id);
                end
            end
        end
        function plotEdge(this, hed0, hed1, color)
            v0 = this.m_mesh.m_heVertexArr(hed0);
            v1 = this.m_mesh.m_heVertexArr(hed1);

            v0_coords = this.m_geom.inputVertexPosition(v0);
            v1_coords = this.m_geom.inputVertexPosition(v1);

            line([v0_coords(1),v1_coords(1)],[v0_coords(2),v1_coords(2)],...
                [v0_coords(3),v1_coords(3)],'linewidth',5,'color',color);
        end

    end
end