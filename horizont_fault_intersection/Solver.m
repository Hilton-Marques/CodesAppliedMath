classdef Solver < Drawer
    properties (Constant)
        m_thershold = 1.5;
    end
    properties            
        m_horizons
        m_faults                
        m_horizon_bb
    end
    methods
        function this = Solver(horizons, faults)
            this = this@Drawer();
            
            this.m_horizons = horizons;
            this.m_faults = faults;

            this.m_horizons.showMesh(this.m_blue);
            this.m_faults.showMesh(this.m_red);
            
            %this.m_horizons.showBoundaryLoops();
            %this.m_faults.showBoundaryLoops();
            this.m_horizon_bb = this.getBB(this.m_horizons.m_geom.m_vertices);
            this.setBB(this.m_horizon_bb);
            this.intersect();
        end
        function intersect(this)
            n_loops = this.m_horizons.m_mesh.m_nboundary_loops;  
            boundaries_projections = {};
            h = [];
            if (n_loops > 0)
                %the boundary loops are at the end of faces array
                first_loop = this.m_horizons.m_mesh.m_nfaces + 1;
                for id_loop = first_loop + 5:(first_loop + n_loops - 1)
                    first_hed = this.m_horizons.m_mesh.m_fHalfEdgeArr(id_loop);
                    prev_hed = first_hed;
                    next_hed = ManifoldSurfaceMesh.INVALID_IND; 
                    search_faces = 1:this.m_faults.m_mesh.m_nfaces;
                    count = 0;
   
                    %this.focusCamera(this.m_horizons.m_boundaries{id_loop - this.m_horizons.m_mesh.m_nfaces});
                    this.focusCamera(this.m_faults, 1:this.m_faults.m_mesh.m_nvertices);

                    view(1,43);
                    view(84,52);
                    boundary_projection = [];                   
                    while (next_hed ~= first_hed)
                        next_hed = this.m_horizons.m_mesh.m_heNextArr(prev_hed);
                        v0 = this.m_horizons.m_mesh.m_heVertexArr(prev_hed);
                        coords_v0 = this.m_horizons.m_geom.inputVertexPosition(v0);
                        plot3(coords_v0(1),coords_v0(2),coords_v0(3),'o','markersize',5,'markerfacecolor','red');
                        %find the ray source
                        next_1 = this.m_horizons.m_mesh.m_heNextArr(this.m_horizons.m_mesh.getTwin(prev_hed));
                        next_2 = this.m_horizons.m_mesh.m_heNextArr(next_1);
                        o = this.m_horizons.m_mesh.m_heVertexArr(next_2);
                        coords_o = this.m_horizons.m_geom.inputVertexPosition(o);
                        %plot3(coords_o(1),coords_o(2),coords_o(3),'o','markersize',5,'markerfacecolor','red');
                        %this.focusCamera([v0,o]);
                        d = coords_v0 - coords_o;
                        %d = d/norm(d);                        
                        [is_intersecting, lambda, triangle_id] = this.m_faults.rayIntersection(coords_o,d,search_faces);
                        h = [h,this.drawRay(coords_o,d)];
                        if is_intersecting
                            k = lambda - 1;
                            p = coords_o + lambda*d;
                            plot3(p(1),p(2),p(3),'o','markersize',5,'markerfacecolor','green');
                            this.m_faults.closestPointToFace(triangle_id, p);
                            this.m_faults.showMesh('cyan',triangle_id);

                            %h = [h,plot3(p(1),p(2),p(3),'o','markersize',5,'markerfacecolor','green')];
                        else
                           delete(h);
                        end
                        pause(0.01);
                        if (is_intersecting && lambda < this.m_thershold)
                            p = coords_o + lambda*d;
                            if count > 0                                 
                                if (triangle_id ~= tri_id_prev)


                                    
                                    %find intermediary point
                                    [isJoinedByEdge, inc] = ...
                                        this.m_faults.m_mesh.isJoinedByEdge(triangle_id, tri_id_prev);
                                    if (isJoinedByEdge)
                                        %approximate rotation
                                        p_k = this.m_faults.project(tri_id_prev,p);
                                        edge = this.m_faults.m_geom.inputVertexPosition(inc);
                                        %plot3(edge(1,1),edge(1,2),edge(1,3),'o','markersize',5,'markerfacecolor','green');
                                        %plot3(edge(2,1),edge(2,2),edge(2,3),'o','markersize',5,'markerfacecolor','green');
                                        p_i = cgeom.segmentSegmentIntersection([p_k;p_prev],edge);
                                    else
                                        %should be joined by vertex
                                        [isJoinedByVertex,v_id] = ...
                                            this.m_faults.m_mesh.isJoinedByVertex(triangle_id, tri_id_prev);
                                        assert(isJoinedByVertex == true);
                                        p_i = this.m_fault.m_geom.inputVertexPosition(v_id);
                                    end
                                    %calculate intermediary point
%                                     p_k = this.m_faults.project(tri_id_prev,p);                                    
%                                     p_k = this.m_faults.closestPointToFace(tri_id_prev, p_k);
                                    %p_i = this.m_faults.closestPointToFace(tri_id_prev, p);
                                    plot3(p_i(1),p_i(2),p_i(3),'o','markersize',5,'markerfacecolor','green');
                                    boundary_projection = [boundary_projection ; p_i];
                                end
                            end     
                            search_faces = this.m_faults.m_mesh.get_adjacentFaces(triangle_id,this.m_faults);
                            %this.m_faults.showMesh('cyan',search_faces);
                            boundary_projection = [boundary_projection ; p];
                            p_prev = p;
                            tri_id_prev = triangle_id;
                            count = count + 1;
                        end
                        prev_hed = next_hed;
                    end
                    boundaries_projections{end+1} = boundary_projection;
                end
                %%show projections
                %plot

            end
        end
        

        function focusCamera(this,obj,  vertices_ids)
            coords = obj.m_geom.inputVertexPosition(vertices_ids);
            bb = this.getBB(coords);
            this.walkTrough(bb);
            this.setBB(bb);

            %this.walkTrough(this.m_parent_bb);
        end
    end
end