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
        function TestNormal(this)
            vertex_id = 1;
            n = this.m_horizon.getVertexNormal(vertex_id);
            hold on
            star_face = this.m_horizon.m_mesh.getStar(vertex_id);
            p0 = this.m_horizon.m_geom.inputVertexPosition(vertex_id);
            this.m_horizon.showMesh(star_face);
            plot3(p0(1), p0(2), p0(3), 'o', 'markerfacecolor','red');
            view(30, 30);
            quiver3(p0(1), p0(2), p0(3), n(1), n(2), n(3), 'color', 'black','linewidth',2);
        end
        function TestExtension(this)
            fac = 1000;
            boundary_loop_id = 1;            
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
                curr_vertex_extended = this.m_horizon.extendVertex(curr_vertex, fac);
                curr_vertex = this.m_horizon.m_geom.inputVertexPosition(curr_vertex);
                tri_1 = [prev_vertex; curr_vertex_extended; prev_vertex_extended];
                tri_2 = [prev_vertex; curr_vertex; curr_vertex_extended];
                trisurf([1,2,3],tri_1(:,1), tri_1(:,2), tri_1(:,3),'FaceColor',"#fdfd96");
                trisurf([1,2,3],tri_2(:,1), tri_2(:,2), tri_2(:,3),'FaceColor',"#fdfd96");
                
                edge_normal_curr = this.m_horizon.getHedNormal(next_hed);
                normal_2 = (edge_normal_prev + edge_normal_curr);
                normal_2 = normal_2 / norm(normal_2);
                edge_curr_extended = curr_vertex + normal_2*fac;

                tri_1 = [prev_vertex; edge_curr_extended; edge_prev_extended];
                tri_2 = [prev_vertex; curr_vertex; edge_curr_extended];
                trisurf([1,2,3],tri_1(:,1), tri_1(:,2), tri_1(:,3),'FaceColor',"#77dd77");
                trisurf([1,2,3],tri_2(:,1), tri_2(:,2), tri_2(:,3),'FaceColor',"#77dd77");
                
                edge_prev_extended = edge_curr_extended;
                edge_normal_prev = edge_normal_curr;


                %plot3(new_vertex(1), new_vertex(2), new_vertex(3), 'o', 'markerfacecolor','red');

                % plot edge and vertices
                prev_hed = next_hed;
                prev_vertex_extended = curr_vertex_extended;
                prev_vertex = curr_vertex;
            end
            a = 1;
        end
    end
    methods (Static)

    end
end