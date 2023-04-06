classdef DFS < handle
    properties
        m_time = 0
        m_prev
        m_pos
        m_mask
        m_list
        m_parent
        m_p;
    end
    methods
        function this = DFS(list)
            figure
            hold on
            n = length(list);
            this.m_list = list;
            this.Plot();            
            this.m_prev = zeros(n,1);
            this.m_pos = zeros(n,1);
            this.m_mask = false(n,1);
            this.m_parent = zeros(n,1);
            for i = 1:n
                if ~this.m_mask(i)
                    this.Visit(i);                    
                end
            end
        end
        function Visit(this, u)
            this.m_time = this.m_time + 1;
            this.m_prev(u) = this.m_time;
            this.m_mask(u) = true;
            highlight(this.m_p, u, 'NodeColor', 'red');
            this.Label(u,this.m_time,'bottom');
            edges = this.m_list{u};
            for v = edges
                if ~this.m_mask(v)
                    this.LabelEdge(u,v,'Tree','black');
                    this.m_parent(v) = u;
                    this.Visit(v);
                else
                    if this.m_prev(u) > this.m_prev(v) && this.m_pos(v) == 0
                        cycle = [u,v];
                        k = this.m_parent(u); 
                        cycle = [cycle,k];
                        while k ~= v
                            k = this.m_parent(k);
                            cycle = [cycle,k];
                        end                        
                        this.DrawCycle(cycle);
                        this.LabelEdge(u,v,'Back','magenta');
                    else
                        this.LabelEdge(u,v,'Forward','magenta');
                    end
                end
            end
            this.m_time = this.m_time + 1;
            this.m_pos(u) = this.m_time;
            this.Label(u,this.m_time,'top');
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