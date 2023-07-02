classdef Fault < handle
    properties
        m_vertices
        m_faces
        m_mesh
        m_geom
        m_mesh_geodesic
        m_algorithm_geodesic
    end
    methods
        function this = Fault(fault_vertices, fault_faces)
            this.m_vertices = fault_vertices;
            this.m_faces = fault_faces;
            this.m_mesh = ManifoldSurfaceMesh(fault_faces);
            this.m_geom = Geometry(fault_vertices);
        end
        function showMesh(this, color, faces_ids)
            if nargin == 2
                faces_ids = 1:this.m_mesh.m_nfaces;
            end
            k = this.m_faces(faces_ids,:);
%             vertices_ids = unique(k(:));
%             vertices = this.m_geom.inputVertexPosition(vertices_ids);
%             trisurf(k, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
%                 'FaceAlpha',0.7,'EdgeColor','none','FaceColor', color);
             trisurf(k, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
                'FaceAlpha',0.7,'FaceColor', color,'EdgeAlpha',0.0);
        end
        function initGeodesicPath(this)
             this.m_mesh_geodesic = geodesic_new_mesh(this.m_vertices, this.m_faces);
             this.m_algorithm_geodesic = geodesic_new_algorithm(this.m_mesh_geodesic, 'exact'); 
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
        function p = project(this, face_id, edge)
            first_hed = this.m_mesh.m_fHalfEdgeArr(face_id);
            v0 = this.m_mesh.m_heVertexArr(first_hed);
            next_hed = this.m_mesh.m_heNextArr(first_hed);
            v1 = this.m_mesh.m_heVertexArr(next_hed);
            next_hed = this.m_mesh.m_heNextArr(next_hed);
            v2 = this.m_mesh.m_heVertexArr(next_hed);
            tri_coords = this.m_geom.inputVertexPosition([v0,v1,v2]);
            p = cgeom.project(tri_coords, edge);
        end
        function p = closestPointToFace(this, face_id, P)
            first_hed = this.m_mesh.m_fHalfEdgeArr(face_id);
            v0 = this.m_mesh.m_heVertexArr(first_hed);
            next_hed = this.m_mesh.m_heNextArr(first_hed);
            v1 = this.m_mesh.m_heVertexArr(next_hed);
            next_hed = this.m_mesh.m_heNextArr(next_hed);
            v2 = this.m_mesh.m_heVertexArr(next_hed);
            tri_coords = this.m_geom.inputVertexPosition([v0,v1,v2]);
            p = cgeom.closestPointToTriangle(tri_coords, P);
        end
        function [is_intersecting, lambda, triangle_id] = rayIntersection(this,o,d,search_faces)
            %for each face, check ray-triangle intersection
            is_intersecting = false;            
            triangle_id = ManifoldSurfaceMesh.INVALID_IND;
            lambda = realmax;
            for i = search_faces
                first_hed = this.m_mesh.m_fHalfEdgeArr(i);
                v0 = this.m_mesh.m_heVertexArr(first_hed);
                next_hed = this.m_mesh.m_heNextArr(first_hed);
                v1 = this.m_mesh.m_heVertexArr(next_hed);
                next_hed = this.m_mesh.m_heNextArr(next_hed);
                v2 = this.m_mesh.m_heVertexArr(next_hed);
                tri_coords = this.m_geom.inputVertexPosition([v0,v1,v2]);
                new_lambda = cgeom.rayTriangleIntersection(o,d,tri_coords);
                if new_lambda < lambda
                    lambda = new_lambda;
                    triangle_id = i;
                end
            end
            if lambda < realmax
                is_intersecting = true;
            end
        end

    end
end