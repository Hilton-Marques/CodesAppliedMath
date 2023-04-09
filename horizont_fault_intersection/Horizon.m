classdef Horizon < handle
    properties
        m_vertices
        m_faces
        m_mesh
        m_geom
        m_boundaries
    end
    methods
        function this = Horizon(horizon_vertices, horizon_faces)
            this.m_vertices = horizon_vertices;
            this.m_faces = horizon_faces;
            this.m_mesh = ManifoldSurfaceMesh(horizon_faces);
            this.m_geom = Geometry(horizon_vertices);
            this.initBoundaryLoops();
        end
        function h = showMesh(this, color, faces_ids)
            if nargin == 2
                faces_ids = 1:this.m_mesh.m_nfaces;
            end
            k = this.m_faces(faces_ids,:);
%             vertices_ids = unique(k(:));
%             vertices = this.m_geom.inputVertexPosition(vertices_ids);
%             trisurf(k, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
%                 'FaceAlpha',0.7,'EdgeColor','none','FaceColor', color);
             h = trisurf(k, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
                'FaceAlpha',0.7,'FaceColor', color,'EdgeAlpha',0.0);
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

        function initBoundaryLoops(this)
            n_loops = this.m_mesh.m_nboundary_loops;
            colors = rand(n_loops,3);

            if (n_loops > 0)
                boundaries = cell(n_loops,1);
                %the boundary loops are at the end of faces array
                first_loop = this.m_mesh.m_nfaces + 1;
                count = 1;
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
                        prev_hed = next_hed;
                    end
                    boundaries{count} = vertices_id;
                    count = count + 1;                    
                end
                this.m_boundaries = boundaries;

            end
        end
        function plotEdge(this, hed0, hed1, color)
            v0 = this.m_mesh.m_heVertexArr(hed0);
            v1 = this.m_mesh.m_heVertexArr(hed1);

            v0_coords = this.m_geom.inputVertexPosition(v0);
            v1_coords = this.m_geom.inputVertexPosition(v1);

            line([v0_coords(1),v1_coords(1)],[v0_coords(2),v1_coords(2)],...
                [v0_coords(3),v1_coords(3)],'linewidth',3,'color',color);
        end

    end
end