classdef ManifPatches < handle
    properties
        m_time = 0
        m_prev
        m_pos
        m_mask_vertices
        m_G
        m_parent
        m_p;
        m_vertices;
        m_faces;
        m_edges;
        m_halfedges;
        m_meshes Mesh;
        m_count_cc = 0;
        m_poses
    end
    methods
        function this = ManifPatches(G, poses)
%             figure
             hold on
%             this.m_vertices = vertices;
%             this.m_faces = faces;
%             this.m_edges = edges;
%             this.m_halfedges = halfedges;
%           
            n = size(G,1);
            this.m_G = G;
            this.m_poses = poses;
            this.m_prev = zeros(n,1);
            this.m_pos = zeros(n,1);
            this.m_mask_vertices = false(n,1);
            
            
            this.m_parent = zeros(n,1);
            for i = 1:n
                if ~this.m_mask_vertices(i)
                    %m = Mesh();  
                    this.m_count_cc = this.m_count_cc + 1;
                    %this.m_meshes(end+1) = m;
                    this.Visit(i,rand(1,3));                    
                end
            end
        end
        function Visit(this, u,color)
            %mesh = this.m_meshes(this.m_count_cc);
            %vertex = this.m_vertices(u);
            %mesh.addVertex(vertex);
            
            this.m_time = this.m_time + 1;
            this.m_prev(u) = this.m_time;
            this.m_mask_vertices(u) = true;

            plot3(this.m_poses(u,1),this.m_poses(u,2),this.m_poses(u,3),'o','markerfacecolor',color,'markersize',10);
            
            
            edges = find(this.m_G(u,:));
            for v = edges
% Mesh update
%                 edge = this.m_edges(v);                 
%                 if ~this.m_edges(edge.id)
%                     mesh.addEdge(edge);
%                     hed_1 = edge.getHed();
%                     mesh.addHalfEdge(hed_1);
%                     if (hed.hasTwin())
%                         hed_2 = edge.getTwin();
%                         mesh.addHalfEdge(hed_2);
%                     end
%                     face = hed_1.getFace();
%                     if ~this.m_faces(face.id)
%                         mesh.addFace(face);
%                         this.m_faces(face.id) = true;
%                     end
%                     this.m_edges(edge.id) = true;
%                 end                
                if ~this.m_mask_vertices(v)
                    this.m_parent(v) = u;
                    this.Visit(v,color);
                end
            end
            this.m_time = this.m_time + 1;
            this.m_pos(u) = this.m_time;

        end
        function Plot(this)
            n = length(this.m_list);
            G = zeros(n,n);
            for i = 1:n
                for j = this.m_list{i}
                    G(i,j) = 1.0;
                end
            end
            if issymmetric(G)
                this.m_p = plot(graph(G),'NodeColor','green','NodeLabel',[],'LineWidth',3.0);
            else
                this.m_p = plot(digraph(G),'NodeColor','green','NodeLabel',[],'LineWidth',1.0);
            end
        end
        function LabelEdge(this,u,v,label,color)
            xi = this.m_p.XData(u);
            yi = this.m_p.YData(u);
            xj = this.m_p.XData(v);
            yj = this.m_p.YData(v);
            x = (xi + xj)/2;
            y = (yi + yj)/2;
            text(x,y,label,'HorizontalAlignment','left','VerticalAlignment','top');
            highlight(this.m_p, [u,v], 'EdgeColor', color);
        end
        function DrawCycle(this,cic)
            x = mean(this.m_p.XData(cic));
            y = mean(this.m_p.YData(cic));
            plot(x,y,'o','markersize',10);
        end
        function Label(this,i,id,format)
            x = this.m_p.XData(i);
            y = this.m_p.YData(i);
            text(x,y,num2str(id),'HorizontalAlignment','left','VerticalAlignment',format);
        end
    end
end