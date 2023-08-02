classdef Tests < handle
    properties
        m_horizon
        m_blue = [0.2588 0.5216 0.9569];
    end
    methods
        function this = Tests()
            file = ReadFile("../meshes/EMB_COMPLETO_160_145_ts.vtk");
            [vertices_horizon, faces_horizon] = file.getMesh;
            this.m_horizon = Horizon(vertices_horizon, faces_horizon);
        end
        function TestNormal(this,vertex_id)
            n = this.m_horizon.getVertexNormal(vertex_id);
            hold on
            star_face = this.m_horizon.m_mesh.getStar(vertex_id);
            p0 = this.m_horizon.m_geom.inputVertexPosition(vertex_id);
            this.m_horizon.showMesh(star_face);
            plot3(p0(1), p0(2), p0(3), 'o', 'markerfacecolor','red');
            view(30, 30);
            quiver3(p0(1), p0(2), p0(3), n(1), n(2), n(3), 'color', 'black','linewidth',2);
        end
        function TestBoundaryNormal(this,vertex_id)
            n = this.m_horizon.getVertexBoundaryNormal(vertex_id);
            hold on
            star_face = this.m_horizon.m_mesh.getStar(vertex_id);
            p0 = this.m_horizon.m_geom.inputVertexPosition(vertex_id);
            this.m_horizon.showMesh(star_face);
            plot3(p0(1), p0(2), p0(3), 'o', 'markerfacecolor','red');
            view(30, 30);
            quiver3(p0(1), p0(2), p0(3), n(1), n(2), n(3), 'color', 'black','linewidth',2);
        end
        
        function TestExtension(this, boundary_loop_id)
            fac = 1000;
            boundary_i = this.m_horizon.m_mesh.m_nfaces + boundary_loop_id;
            faces_id = this.m_horizon.m_mesh.getFacesFromBoundary(boundary_i);
            hold on
            view(30,30);
            h = this.m_horizon.showMesh(faces_id, this.m_blue);
            
            % extend vertices
            first_hed = this.m_horizon.m_mesh.m_fHalfEdgeArr(boundary_i);
            prev_hed = first_hed;

            prev_vertex = this.m_horizon.m_mesh.m_heVertexArr(prev_hed);
            prev_vertex_extended = this.m_horizon.extendVertex(prev_vertex, fac);
            prev_vertex = this.m_horizon.m_geom.inputVertexPosition(prev_vertex);
            curr_vertex = prev_vertex;

            edge_normal_prev = this.m_horizon.getHedNormal(prev_hed);
            edge_prev_extended = prev_vertex_extended;
            
            next_hed = ManifoldSurfaceMesh.INVALID_IND;            
            vertices_id = [];            
            while (next_hed ~= first_hed)
                next_hed = this.m_horizon.m_mesh.m_heNextArr(prev_hed);
                curr_vertex = this.m_horizon.m_mesh.m_heVertexArr(next_hed);

                vertices_id = [vertices_id, curr_vertex];
                
                %create two new triangles   

                %method one with cotan formula
                curr_vertex_extended = this.m_horizon.extendVertex(curr_vertex, fac);
                curr_vertex = this.m_horizon.m_geom.inputVertexPosition(curr_vertex);
                tri_1 = [prev_vertex; curr_vertex_extended; prev_vertex_extended];
                tri_2 = [prev_vertex; curr_vertex; curr_vertex_extended];
                %trisurf([1,2,3],tri_1(:,1), tri_1(:,2), tri_1(:,3),'FaceColor',"#fdfd96");
                %trisurf([1,2,3],tri_2(:,1), tri_2(:,2), tri_2(:,3),'FaceColor',"#fdfd96");
                
                %method two with mean edge normal formula
                edge_normal_curr = this.m_horizon.getHedNormal(next_hed);
                normal_2 = (edge_normal_prev + edge_normal_curr);
                normal_2 = normal_2 / norm(normal_2);
                edge_curr_extended = curr_vertex + normal_2*fac;
                

                tri_1 = [prev_vertex; edge_curr_extended; edge_prev_extended];
                tri_2 = [prev_vertex; curr_vertex; edge_curr_extended];
                
                trisurf([1,2,3],tri_1(:,1), tri_1(:,2), tri_1(:,3),'FaceColor',"#77dd77");
                trisurf([1,2,3],tri_2(:,1), tri_2(:,2), tri_2(:,3),'FaceColor',"#77dd77");
                
                %plot3(new_vertex(1), new_vertex(2), new_vertex(3), 'o', 'markerfacecolor','red');

                % plot edge and vertices
                prev_hed = next_hed;

                % method one updates
                prev_vertex_extended = curr_vertex_extended;
                prev_vertex = curr_vertex;

                % method two updates
                edge_prev_extended = edge_curr_extended;
                edge_normal_prev = edge_normal_curr;
            end
        end
        function TestRetriangulation(this)
            file = ReadFile("../meshes/result_triangles_6.vtk");
            [V, F] = file.getMesh;
            hold on;
            trisurf(F,V(:,1),V(:,2),V(:,2));
            trisurf(F,V(:,1),V(:,2),0*V(:,2));
            horizon = Horizon(V,F);
            bnd = horizon.m_boundaries{1}';
            n = size(bnd,1);    
            shiftedVertices = circshift(bnd, -1);
            edges = [bnd, shiftedVertices];
            DT = delaunayTriangulation(V(:,1),V(:,2),double(edges));
            SF = DT.ConnectivityList;
            [SV, SF] = triangle(V(:,1:2), edges, []);
            trisurf(SF,V(:,1),V(:,2),V(:,2));
            mesh = {struct('vertices', V, 'faces',SF)};
            horizon.exportTs(mesh, 'delaunau_fault_six.vtk');
            a = 1;
        end
        function TestSDFExtension(this, boundary_loop_id)
            [V,F] = this.m_horizon.createPatch(boundary_loop_id);
            %this.objWriter(V,F,"patch.obj");
            %[V,F] = load_mesh('model.obj');
           
            % Extract offset at minus 3% of bounind box diagonal length
            iso = 0.03;
            % Resolution grid → resolution of extracted offset surface
            side = 60;
            % Amount of smoothing to apply to distance field
            sigma = 1.4;
            bbd = norm(max(V)-min(V));
            % Pad generously for positive iso values
%             [BC,side,r] = voxel_grid([V;max(V)+iso*1;min(V)-iso*1],side);
%             D = signed_distance(BC,V,F,'SignedDistanceType','unsigned');
%             D = reshape(D,side([2 1 3]));
%             % Smooth signed distance field
%             D = imfilter(D,fspecial('gaussian',9,sigma),'replicate');
%             BC3 = reshape(BC,[side([2 1 3]) 3]);
%             % Use matlab's built-in marching cubes iso-surface mesher (works most of the time)
%             surf = isosurface(BC3(:,:,:,1),BC3(:,:,:,2),BC3(:,:,:,3),D,100);
%             SV = surf.vertices;
%             SF = surf.faces;
            [SV, SF] = signed_distance_isosurface(V,F,'Level',1000,'SignedDistanceType','unsigned');
            clf;
            hold on;
            t  = tsurf(F,V,'EdgeColor','none',fsoft,  'FaceVertexCData',repmat(blue,size(V,1),1),'FaceAlpha',1+(iso<0)*(0.35-1),fphong);
            ts = tsurf(SF,SV,'EdgeAlpha',0.2+(iso<0)*(0-0.2),fsoft,'FaceVertexCData',repmat(orange,size(SV,1),1),fphong,'FaceAlpha',1+(iso>0)*(0.2-1));
            apply_ambient_occlusion(ts);
            hold off;
            axis equal;
            view(-20,20)
            camlight;
            t.SpecularStrength = 0.04;
            l = light('Position',[5 -5 10],'Style','local');
            %add_shadow(t,l,'Color',0.8*[1 1 1],'Fade','local','Ground',[0 0 -1 min([V(:,3);SV(:,3)])]);
            set(gca,'pos',[0 0 1 1])
            set(gca,'Visible','off');
            set(gcf,'Color','w');
            drawnow;
        end
        function FF = RemoveBadTriangle(this, eps)
            file = ReadFile("../meshes/result_triangles.vtk");
            [V, F] = file.getMesh;
            [FF] = collapse_small_triangles(V, F, eps);
%             verts_ids = sort(unique(FF(:)));
%             VV = V(verts_ids, :);
%             vertice_mapping = zeros(max(verts_ids),1);
%             vertice_mapping(verts_ids) = 1:size(verts_ids,1);
%             FF = vertice_mapping(FF);
%             mesh = {struct('vertices', VV, 'faces',FF)};
%             Horizon.exportTs(mesh, "new_result_2.vtk");
        end
        function FF = MyRemoveBadTriangle(this, eps)
            file = ReadFile("../meshes/result_triangles_0108.vtk");
            [V, F] = file.getMesh;
            [FF] = cgeom.collapse_small_areas(V,F,eps);
            verts_ids = sort(unique(FF(:)));
            VV = V(verts_ids, :);
            vertice_mapping = zeros(max(verts_ids),1);
            vertice_mapping(verts_ids) = 1:size(verts_ids,1);
            FF = vertice_mapping(FF);
            mesh = {struct('vertices', VV, 'faces',FF)};
            Horizon.exportTs(mesh, "new_result_2.vtk");
        end
    end
    methods (Static)
        function objWriter(vertices,faces, filename)
            new_filename = split(filename, '.') ;
            new_filename = strcat(new_filename{1} , '.obj');

            % Open the output file for writing
            fileID = fopen(new_filename, 'w');

            % Write the vertex data to the file
            fprintf(fileID, 'v %f %f %f\n', vertices');

            % Write the face data to the file
            fprintf(fileID, 'f %d %d %d\n', faces');
            fclose(fileID);
        end
    end
end