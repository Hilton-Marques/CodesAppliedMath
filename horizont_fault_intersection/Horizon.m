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
        function new_vertex = extendVertex(this, id, fac)
            normal = this.getVertexNormal(id);
            pi = this.m_geom.inputVertexPosition(id);
            new_vertex = pi + fac * normal;
        end
        function h = showMesh(this, faces_ids, color)            
            if nargin == 1
                faces_ids = 1:this.m_mesh.m_nfaces;
                color = "#A7C7E7";
            end
            if nargin == 2
                color = "#A7C7E7";
            end
            k = this.m_faces(faces_ids,:);
%             vertices_ids = unique(k(:));
%             vertices = this.m_geom.inputVertexPosition(vertices_ids);
%             trisurf(k, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
%                 'FaceAlpha',0.7,'EdgeColor','none','FaceColor', color);
             h = trisurf(k, this.m_vertices(:,1), this.m_vertices(:,2), this.m_vertices(:,3),...
                'FaceAlpha',0.7,'FaceColor', color,'EdgeAlpha',1.0);
        end
        function normal = getVertexNormal(this, id)
            star_face = this.m_mesh.getStar(id);
            vertex_ids = this.m_faces(star_face,:);
            p0 = this.m_geom.inputVertexPosition(id);

            n_triangles = size(vertex_ids,1);
            normal = zeros(1,3);            
            for i = 1:n_triangles                
                for id_face = 1:3
                    face_id = vertex_ids(i, id_face);
                    if (face_id == id)
                        break;
                    end
                end                
                p1 = this.m_geom.inputVertexPosition(vertex_ids(i, mod(id_face, 3) + 1));
                p2 = this.m_geom.inputVertexPosition(vertex_ids(i, mod(id_face + 1, 3) + 1));
                %plot3(p1(1), p1(2), p1(3), 'o', 'markerfacecolor','#77dd77');
                %plot3(p2(1), p2(2), p2(3), 'o', 'markerfacecolor','#77dd77');
                %m = (p0 + p1 + p2) / 3;
                l1 = p1 - p0;
                l2 = p2 - p0;
                l3 = p1 - p2;
                n = cross(l2, l1);
                n = n / norm(n);
                %quiver3(m(1), m(2), m(3), n(1), n(2), n(3), 'color', '#77dd77');
                grad = 0.5 * cross(n , l3);
                %quiver3(p0(1), p0(2), p0(3), grad(1), grad(2), grad(3), 'color', '#ff6961');
                normal = normal + grad;
            end
            normal = normal/norm(normal);
        end
        function edge_normal = getHedNormal(this, hed_id)
            v0 = this.m_mesh.m_heVertexArr(hed_id);
            v1 = this.m_mesh.m_heVertexArr(this.m_mesh.m_heNextArr(hed_id));
            edge = this.m_geom.inputVertexPosition(v1) - this.m_geom.inputVertexPosition(v0);
            face_id = this.m_mesh.m_heFaceArr(this.m_mesh.getTwin(hed_id));
            vertices = this.m_faces(face_id,:);
            pts = this.m_geom.inputVertexPosition(vertices);
            n = cross(pts(2,:) - pts(1,:), pts(3,:) - pts(1,:));
            edge_normal = cross(n, edge);
            edge_normal = edge_normal/norm(edge_normal);
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