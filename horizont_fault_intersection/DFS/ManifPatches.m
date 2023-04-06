classdef ManifPatches < handle
    properties
        m_time = 0
        m_prev
        m_pos
        m_mask_vertices
        m_list
        m_parent
        m_p;
        m_vertices;
        m_faces;
        m_edges;
        m_halfedges;
        m_meshes Mesh;
        m_count_cc = 0;
    end
    methods
        function this = ManifPatches(list, vertices, faces, edges, halfedges)
            figure
            hold on
            n = length(list);
            this.m_list = list;
            this.Plot();            
            this.m_vertices = vertices;
            this.m_faces = faces;
            this.m_edges = edges;
            this.m_halfedges = halfedges;
            this.m_prev = zeros(n,1);
            this.m_pos = zeros(n,1);
            this.m_mask_vertices = false(n,1);
            n_edges = size(edges,2);
            n_faces = size(faces,2);
            this.m_mask_edges = false(n_edges,1);
            this.m_mask_faces = false(n_faces,1);
            this.m_parent = zeros(n,1);
            for i = 1:n
                if ~this.m_mask(i)
                    m = Mesh();  
                    this.m_count_cc = this.m_count_cc + 1;
                    this.m_meshes(end+1) = m;
                    this.Visit(i);                    
                end
            end
        end
        function Visit(this, u)
            mesh = this.m_meshes(this.m_count_cc);
            vertex = this.m_vertices(u);
            mesh.addVertex(vertex);
            
            this.m_time = this.m_time + 1;
            this.m_prev(u) = this.m_time;
            this.m_mask_vertices(u) = true;
            highlight(this.m_p, u, 'NodeColor', 'red');
            this.Label(u,this.m_time,'bottom');
            
            edges = this.m_list{u};
            for v = edges
                edge = this.m_edges(v);                 
                if ~this.m_edges(edge.id)
                    mesh.addEdge(edge);
                    hed_1 = edge.getHed();
                    mesh.addHalfEdge(hed_1);
                    if (hed.hasTwin())
                        hed_2 = edge.getTwin();
                        mesh.addHalfEdge(hed_2);
                    end
                    face = hed_1.getFace();
                    if ~this.m_faces(face.id)
                        mesh.addFace(face);
                        this.m_faces(face.id) = true;
                    end
                    this.m_edges(edge.id) = true;
                end                
                if ~this.m_mask_vertices(v)
                    this.m_parent(v) = u;
                    this.Visit(v);
                end
            end
            this.m_time = this.m_time + 1;
            this.m_pos(u) = this.m_time;
            this.Label(u,this.m_time,'top');
            highlight(this.m_p, u, 'NodeColor', 'black');
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