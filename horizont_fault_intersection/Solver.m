classdef Solver < Drawer
    properties (Constant)
        m_thershold = 1.25;
    end
    properties            
        m_horizons
        m_faults                
        m_horizon_bb
        m_boundary_handles = [];
    end
    methods
        function this = Solver(horizons, faults)
            this = this@Drawer();
            
            this.m_horizons = horizons;
            this.m_faults = faults;
            %this.m_faults.initGeodesicPath();
            
            this.m_faults.showMesh(this.m_red);
            this.poissonReconstruction(this.m_faults);
            for i = 1:this.m_horizons.m_mesh.m_nboundary_loops
                boundary_i = this.m_horizons.m_mesh.m_nfaces + i;
                faces_id = this.m_horizons.m_mesh.getFacesFromBoundary(boundary_i);
                h = this.m_horizons.showMesh(faces_id,this.m_blue);
                this.m_boundary_handles = [this.m_boundary_handles, h];
               % this.Extend(5000, i);
            end
%             this.m_horizon_bb = this.getBB(this.m_horizons.m_geom.m_vertices);
%             this.setBB(this.m_horizon_bb);
%%Take 1
            this.m_horizons.showMesh();
            this.m_horizon_bb = this.getBB(this.m_horizons.m_geom.m_vertices);
            this.setBB(this.m_horizon_bb);
            %ids = unique(data(:));
            %h = this.m_horizons.showMesh(ids,this.m_red);
            %this.focusCamera(this.m_horizons,  ids);
%             for i = 1:90
%                 this.update();
%             end
%             this.save();
            
%% Take 2
            % this.m_horizons.showBoundaryLoops();
%             for i = 1:90
%                 this.update();
%             end
%             this.save('take_2.gif');
            %this.m_faults.showBoundaryLoops();
            this.m_faults.showMesh(this.m_red);
%% Take 3
            %view(-17, 57);
            %this.exportFrame("take_3");
%% Take 4
            %this.showGeodesicPath()
            this.intersect();
        end
        function showNewTriangles()
            triangle_id = 3;
            
        end
        function Extend(this, fac, boundary_loop_id)
            boundary_i = this.m_horizons.m_mesh.m_nfaces + boundary_loop_id;
            % extend vertices
            first_hed = this.m_horizons.m_mesh.m_fHalfEdgeArr(boundary_i);
            prev_hed = first_hed;

            prev_vertex = this.m_horizons.m_mesh.m_heVertexArr(prev_hed);
            prev_vertex_extended = this.m_horizons.extendVertex(prev_vertex, fac);
            prev_vertex = this.m_horizons.m_geom.inputVertexPosition(prev_vertex);
            curr_vertex = prev_vertex;
            
            next_hed = ManifoldSurfaceMesh.INVALID_IND;            
            vertices_id = [];            
            while (next_hed ~= first_hed)
                next_hed = this.m_horizons.m_mesh.m_heNextArr(prev_hed);
                curr_vertex = this.m_horizons.m_mesh.m_heVertexArr(next_hed);

                vertices_id = [vertices_id, curr_vertex];
                
                %create two new triangles   

                %method one with cotan formula
                curr_vertex_extended = this.m_horizons.extendVertex(curr_vertex, fac);
                curr_vertex = this.m_horizons.m_geom.inputVertexPosition(curr_vertex);
                tri_1 = [prev_vertex; curr_vertex_extended; prev_vertex_extended];
                tri_2 = [prev_vertex; curr_vertex; curr_vertex_extended];
                trisurf([1,2,3],tri_1(:,1), tri_1(:,2), tri_1(:,3),'FaceColor',"#fdfd96");
                trisurf([1,2,3],tri_2(:,1), tri_2(:,2), tri_2(:,3),'FaceColor',"#fdfd96");
                
        
                prev_hed = next_hed;

                prev_vertex_extended = curr_vertex_extended;
                prev_vertex = curr_vertex;

            end
            a = 1;
        end
        function intersect(this)
            n_loops = this.m_horizons.m_mesh.m_nboundary_loops;  
            boundaries_projections = {};
            h = [];
            if (n_loops > 0)
                %the boundary loops are at the end of faces array
                first_loop = this.m_horizons.m_mesh.m_nfaces + 1;
                for id_loop = first_loop + 1:(first_loop + n_loops - 1)
                    first_hed = this.m_horizons.m_mesh.m_fHalfEdgeArr(id_loop);
                    prev_hed = first_hed;
                    next_hed = ManifoldSurfaceMesh.INVALID_IND; 
                    search_faces = 1:this.m_faults.m_mesh.m_nfaces;
                    count = 0;
                    initialized = false;
                    %this.focusCamera(this.m_horizons.m_boundaries{id_loop - this.m_horizons.m_mesh.m_nfaces});
                    %Take 6
%                     this.focusCamera(this.m_faults, 1:this.m_faults.m_mesh.m_nvertices);
%                     for i = 1:n_loops
%                         if (i == 6)
%                             continue;
%                         end
%                         delete(this.m_boundary_handles(i));
%                     end
%                     view(1,43);
%                     view(84,52);
                    boundary_projection = [];     
                    lambdas = [];
                    take_5 = true;
                    while (next_hed ~= first_hed)
                        next_hed = this.m_horizons.m_mesh.m_heNextArr(prev_hed);
                        v0 = this.m_horizons.m_mesh.m_heVertexArr(prev_hed);
                        coords_v0 = this.m_horizons.m_geom.inputVertexPosition(v0);
                        %plot3(coords_v0(1),coords_v0(2),coords_v0(3),'o','markersize',3,'markerfacecolor','black');
                        %find the ray source
                        next_1 = this.m_horizons.m_mesh.m_heNextArr(this.m_horizons.m_mesh.getTwin(prev_hed));
                        next_2 = this.m_horizons.m_mesh.m_heNextArr(next_1);
                        o = this.m_horizons.m_mesh.m_heVertexArr(next_2);
                        coords_o = this.m_horizons.m_geom.inputVertexPosition(o);
                        %plot3(coords_o(1),coords_o(2),coords_o(3),'o','markersize',3,'markerfacecolor','yellow');
                        %this.focusCamera([v0,o]);
                        d = coords_v0 - coords_o;
                        %d = d/norm(d);                        
                        [is_intersecting, lambda, triangle_id] = this.m_faults.rayIntersection(coords_o,d,search_faces);
                        if take_5
                            h_i = this.drawRay(coords_o,d);
                            h = [h,h_i];
                            text_pose = coords_o + 3*d;
                            if lambda == realmax
                                h_i = text(text_pose(1),text_pose(2),text_pose(3),...
                                    strcat('$\lambda',{''}, '=',{' '},'\infty$'),...
                                    'interpreter','latex','FontSize',16,'VerticalAlignment','bottom');
                                h = [h,h_i];
                            else
                                h_i = text(text_pose(1),text_pose(2),text_pose(3),...
                                    strcat('$\lambda',{''}, '=',{' '},num2str(lambda,2),'$'),...
                                    'interpreter','latex','FontSize',16,'VerticalAlignment','bottom');
                                h = [h,h_i];
                            end
                            if ~initialized
                                this.update(1);
                            else
                                this.get();
                            end
                            delete(h);
                        end
                        prev_hed = next_hed;
                        if (is_intersecting && lambda < this.m_thershold)
                            if count < 0
                                initialized = true;
                                this.focusCamera(this.m_faults, 1:this.m_faults.m_mesh.m_nvertices);
                                for i = 1:n_loops
                                    if (i == 6)
                                        continue;
                                    end
                                    delete(this.m_boundary_handles(i));
                                end
                                view(84,52);
                            end
                            p = coords_o + lambda*d;
                            plot3(p(:,1), p(:,2), p(:,3),'o','MarkerSize',4,'MarkerFaceColor','green','color','green');                            
                            lambdas = [lambdas, lambda];
                            if count < 0                                 
                                if (triangle_id ~= tri_id_prev)
                                    %create geodesic path
                                    source_point{1} = geodesic_create_surface_point('face',tri_id_prev,p_prev);
                                    geodesic_propagate(this.m_faults.m_algorithm_geodesic, source_point); 
                                    destinations = geodesic_create_surface_point('face',triangle_id,p);
                                    path = geodesic_trace_back(this.m_faults.m_algorithm_geodesic,...
                                        destinations);     %find a shortest path from source to destination
                                    [x,y,z] = extract_coordinates_from_path(path);
                                    p_i = [flip(x(2:end-1)), flip(y(2:end-1)), flip(z(2:end-1))]; %exclude source and dest points
                                    %p_i = flip(p_i);
                                    plot3(p_i(:,1), p_i(:,2), p_i(:,3),'o','MarkerSize',5,'MarkerFaceColor','cyan','color','cyan');
                                    boundary_projection = [boundary_projection ; p_i];
                                end
                            end     
                            %search_faces = this.m_faults.m_mesh.get_adjacentFaces(triangle_id,this.m_faults);
                            %this.m_faults.showMesh('cyan',search_faces);
                            boundary_projection = [boundary_projection ; p];
                            p_prev = p;
                            tri_id_prev = triangle_id;
                            count = count + 1;
                        end
                    end
                    if ~isempty(boundary_projection)
                        %Take 5 
                        this.save('take_5.gif',0.2);
                        %Take 6
                        %this.exportFrame('take_6');
                        %Take 7
                        %this.exportFrame('take 7');
                        %Take 8
                        this.get();
                        n = size(boundary_projection,1);
                        for i = 1:n-1
                            p_i = boundary_projection(i,:);
                            p_j = boundary_projection(i+1,:);
                            line([p_i(1,1), p_j(1,1)],[p_i(1,2), ...
                                p_j(1,2)],[p_i(1,3), p_j(1,3)],'linewidth',2.5,'color','black');
                            this.get();
                            %pause(0.01);
                            %drawnow;
                        end
                        this.save('take_8.gif', 0.1);
                    end
                    boundaries_projections{end+1} = boundary_projection;
                end
                %%show projections
                %plot
            end
        end
        function [vertices,faces] = extendHorizons(this, fac)
            vertices = this.m_horizons.m_vertices;
            faces = this.m_horizons.m_faces;
            n_verts = size(vertices,2);
            for i = 1:this.m_horizons.m_mesh.m_nboundary_loops
                boundary_i = this.m_horizons.m_mesh.m_nfaces + i;
                faces_id = this.m_horizons.m_mesh.getFacesFromBoundary(boundary_i);
                h = this.m_horizons.showMesh(faces_id,this.m_blue);
                this.m_boundary_handles = [this.m_boundary_handles, h];
                [new_vertices, new_faces] = this.Extend(fac, i,n_verts);
                n_verts = size(new_vertices,2) + n_verts;
                vertices = [vertices; new_vertices];
                faces = [faces; new_faces];
            end
        end

      function showGeodesicPath(this)
            axis tight
            source_id = 1;
            source_pt = this.m_faults.m_geom.inputVertexPosition(source_id);
            source_point{1} = geodesic_create_surface_point('vertex',...
                source_id,source_pt);
            geodesic_propagate(this.m_faults.m_algorithm_geodesic, source_point);
            dest_id = this.m_faults.m_mesh.m_nvertices;
            dest_pt = this.m_faults.m_geom.inputVertexPosition(dest_id);
            destinations = geodesic_create_surface_point('vertex',...
                dest_id, dest_pt);
            path = geodesic_trace_back(this.m_faults.m_algorithm_geodesic,...
                destinations);     %find a shortest path from source to destination
            [x,y,z] = extract_coordinates_from_path(path);
            path = flip([x,y,z])
            x = path(:,1); y = path(:,2); z = path(:,3);
            view(86,16);
            this.get();
            plot3(source_pt(1),source_pt(2),source_pt(3),'o','markersize',10,...
                'MarkerFaceColor',this.m_blue);   
            this.get();
            plot3(dest_pt(1),dest_pt(2),dest_pt(3),'o','markersize',10,...
                'MarkerFaceColor','green');
            this.get();
            for i = 1:size(x,1) - 1
                line([x(i), x(i+1)], [y(i), y(i+1)], [z(i), z(i+1)],'LineWidth',2,'color','black'); 
                if (i ~= size(x,1) - 1)
                    plot3(x(i+1),y(i+1),z(i+1),'o','markersize',4,...
                        'MarkerFaceColor','cyan','color','cyan');
                end
                this.get();
            end
            %plot3(x,y,z,'LineWidth',2,'color','black');
            this.save("take_4",0.25);
            %this.exportFrame("take_4");
        end
        function poissonReconstruction(this, manifold_surface)
            points = manifold_surface.m_geom.inputVertexPosition();
            %faces = delaunay(points(:,1),points(:,2), points(:,3));
            %trisurf(faces,points(:,1),points(:,2), points(:,3));
            ptCloud = pointCloud(points);
            %plot3(points(:,1), points(:,2), points(:,3),'o','markerSize',5)
            [mesh,depth,perVertexDensity] = pc2surfacemesh(ptCloud,"poisson",4);
            
            trisurf(mesh.Faces, mesh.Vertices(:,1), mesh.Vertices(:,2), mesh.Vertices(:,3), 'FaceColor', 'cyan');
            %pcshow(ptCloud);
            %surfaceMeshShow(mesh);
            ptCloud = pcread("teapot.ply");
            mesh = pc2surfacemesh(ptCloud,"poisson");
            surfaceMeshShow(mesh)
            a = 1;
        end

        function focusCamera(this,obj,  vertices_ids, withWalkThrough)
            if nargin == 3
                withWalkThrough = false;
            end
            coords = obj.m_geom.inputVertexPosition(vertices_ids);
            bb = this.getBB(coords);
            if withWalkThrough
                this.walkTrough(bb);
            end
            this.setBB(bb);

            %this.walkTrough(this.m_parent_bb);
        end
    end
end