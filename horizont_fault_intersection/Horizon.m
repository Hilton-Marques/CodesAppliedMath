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
            %this.exportPatches([1,5], "patch_1_5.vtk");
            %this.exportPatches([1,5], "patch_1_5.ts");

        end
        function new_vertex = extendVertex(this, id, fac)
            %normal = this.getVertexNormal(id);
            normal = this.getVertexBoundaryNormal(id);
            pi = this.m_geom.inputVertexPosition(id);
            new_vertex = pi + fac * normal;
        end
        function n = getVertexBoundaryNormal(this, vertex_id)
            star_face = this.m_mesh.getStar(vertex_id);
            first_face = star_face(1);
            n_i = this.getFaceNormal(first_face);            
            last_face = star_face(end);
            n_j = this.getFaceNormal(last_face);
            hed = this.m_mesh.m_fHalfEdgeArr(first_face);
            while (true)
                curr_id = this.m_mesh.m_heVertexArr(hed);
                if (curr_id == vertex_id)
                    break
                end
                hed = this.m_mesh.m_heNextArr(hed);
            end
            v_i = this.m_mesh.m_heVertexArr(this.m_mesh.getTwin(hed));
            hed = this.m_mesh.m_fHalfEdgeArr(last_face);
            while (true)
                curr_id = this.m_mesh.m_heVertexArr(this.m_mesh.getTwin(hed));
                if (curr_id == vertex_id)
                    break
                end
                hed = this.m_mesh.m_heNextArr(hed);
            end
            v_j = this.m_mesh.m_heVertexArr(hed);
            %calculate normal
            pi = this.m_geom.inputVertexPosition(v_i);
            p0 = this.m_geom.inputVertexPosition(vertex_id);
            pj = this.m_geom.inputVertexPosition(v_j);
            l1 = cross(n_i, p0 - pi);
            l1 = l1/norm(l1);
            l2 = cross(n_j, pj - p0);
            l2 = l2/norm(l2);
            n = (l1 + l2);
            n = n / norm(n);           
        end
        function n = getFaceNormal(this, face_id)
            vertices = this.m_faces(face_id,:);            
            pts = this.m_geom.inputVertexPosition(vertices);
            n = cross(pts(2,:) - pts(1,:), pts(3,:) - pts(1,:));
            %n = n/norm(n);            
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
            
%             hold on
%             star_face = this.m_mesh.getStar(id);
%             p0 = this.m_geom.inputVertexPosition(id);
%             this.showMesh(star_face);
%             plot3(p0(1), p0(2), p0(3), 'o', 'markerfacecolor','red');
%             view(30, 30);

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
                m = (p0 + p1 + p2) / 3;
                l1 = p1 - p0;
                l2 = p2 - p0;
                l3 = p1 - p2;
                n = cross(l2, l1);
                %n = n / norm(n);
                n = n;
                %quiver3(m(1), m(2), m(3), n(1), n(2), n(3), 'color', '#77dd77');
                grad = 0.5 * cross(n , l3);
                %quiver3(p0(1), p0(2), p0(3), grad(1), grad(2), grad(3), 'color', '#ff6961');
                normal = normal + grad;
            end
            normal = normal/norm(normal);
            %normal = normal;
            %hold off

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
        function exportPatches(this, ids, filename)
            num_meshes = size(ids,2);
            meshes = cell(1,num_meshes);
            n_verts = 0;
            for i = 1:num_meshes
                [vertices, triangles] = this.createPatch(ids(i));
                meshes{i} = struct('vertices', vertices, 'faces',triangles);
            end
            %this.exportVTK(meshes, filename);
            this.exportTs(meshes, filename);
        end
        function [vertices, triangles] = createPatch(this, id)
            boundary_i = this.m_mesh.m_nfaces + id;
            mesh = this.m_mesh.getFacesFromBoundary(boundary_i);
            mesh = this.m_faces(mesh,:);
            verts_ids = sort(unique(mesh(:)));
            vertices = this.m_geom.inputVertexPosition(verts_ids);
            vertice_mapping = zeros(max(verts_ids),1);
            vertice_mapping(verts_ids) = 1:size(verts_ids,1);
            triangles = vertice_mapping(mesh);

            
        end

    end
    methods(Static)
        function exportTs(meshes, filename)
            % Open the file for writing
            fid = fopen(filename, 'w');

            % Write the header
            fprintf(fid, 'Region_1 0\n');
            fprintf(fid, 'TFACE\n');

            % Write the mesh data
            num_meshes = numel(meshes);
            total_vertices = sum(cellfun(@(x) size(x.vertices, 1), meshes));
            total_faces = sum(cellfun(@(x) size(x.faces, 1), meshes));            
            for i = 1:num_meshes
                vertices = meshes{i}.vertices;
                num_vertices = size(vertices, 1);
                fprintf(fid, 'VRTX %f %f %f\n', vertices.');
            end
            fprintf(fid, '\n');
            fprintf(fid, 'CELLS %d %d\n', total_faces, total_faces * 4);
            offset = 0;
            for i = 1:num_meshes
                vertices = meshes{i}.vertices;
                num_vertices = size(vertices, 1);
                faces = meshes{i}.faces;
                num_faces = size(faces, 1);
                cell_data = [repmat(3, num_faces, 1) faces - 1 + offset];
                fprintf(fid, '%d %d %d %d\n', cell_data.');
                offset = offset + num_vertices;
            end
            fprintf(fid, '\n');
            fprintf(fid, 'CELL_TYPES %d\n', total_faces);
            for i = 1:num_meshes
                faces = meshes{i}.faces;
                num_faces = size(faces, 1);
                cell_types = repmat(5, num_faces, 1);
                fprintf(fid, '%d\n', cell_types);
            end

            % Close the file
            fclose(fid);
        end
        function exportVTK(meshes, filename)
            % Open the file for writing
            fid = fopen(filename, 'w');

            % Write the header
            fprintf(fid, '# vtk DataFile Version 2.0\n');
            fprintf(fid, 'Unstructured Grid\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

            % Write the mesh data
            num_meshes = numel(meshes);
            total_vertices = sum(cellfun(@(x) size(x.vertices, 1), meshes));
            total_faces = sum(cellfun(@(x) size(x.faces, 1), meshes));
            fprintf(fid, 'POINTS %d double\n', total_vertices);
            for i = 1:num_meshes
                vertices = meshes{i}.vertices;
                num_vertices = size(vertices, 1);
                fprintf(fid, '%f %f %f\n', vertices.');
            end
            fprintf(fid, '\n');
            fprintf(fid, 'CELLS %d %d\n', total_faces, total_faces * 4);
            offset = 0;
            for i = 1:num_meshes
                vertices = meshes{i}.vertices;
                num_vertices = size(vertices, 1);
                faces = meshes{i}.faces;
                num_faces = size(faces, 1);
                cell_data = [repmat(3, num_faces, 1) faces - 1 + offset];
                fprintf(fid, '%d %d %d %d\n', cell_data.');
                offset = offset + num_vertices;
            end
            fprintf(fid, '\n');
            fprintf(fid, 'CELL_TYPES %d\n', total_faces);
            for i = 1:num_meshes
                faces = meshes{i}.faces;
                num_faces = size(faces, 1);
                cell_types = repmat(5, num_faces, 1);
                fprintf(fid, '%d\n', cell_types);
            end

            % Close the file
            fclose(fid);
        end
        function exportMeshesToTsurf(meshes, outputPath)
            % Open the output file for writing
            fileID = fopen(outputPath, 'w');

            % Write the header
            fprintf(fileID, 'GOCAD TSurf 1.0\n');

            % Loop over the meshes and write them to the file
            for i = 1:numel(meshes)
                % Write the mesh name
                fprintf(fileID, 'NAME %s\n', meshes{i}.name);

                % Write the vertex count
                vertexCount = size(meshes{i}.vertices, 1);
                fprintf(fileID, 'VRTX %d %f %f %f\n', [1:vertexCount; meshes{i}.vertices.']);

                % Write the face count
                faceCount = size(meshes{i}.faces, 1);
                fprintf(fileID, 'TRGL %d %d %d\n', meshes{i}.faces.');
            end

            % Close the file
            fclose(fileID);
        end
    end
end