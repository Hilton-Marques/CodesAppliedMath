classdef SparseMatrix < handle
    properties
        m_M;
        m_G
        m_size
        m_P
        m_E;
    end
    methods
        function this = SparseMatrix(M)
            this.m_M = M;
            this.m_G = digraph(M);
            this.m_size = size(M);
            %init pose
            figure
            axis off
            hold on;
            M = M - diag(diag(M));
            if (M == M')
                G = graph(M);
            else
            G = digraph(M);
            end
            g = this.plot(G);
            this.m_P = [g.XData; g.YData];
            this.m_E = eye(this.m_size);
        end
        %plot without diagonal self loops
        function g = show(this, T)
            if (nargin == 1)
                T = this.m_M;
            end
            M = T - diag(diag(T));
            G = digraph(M);
            g = this.plot(G);
            this.SetGraphPose(g);
        end
        function RemoveNode(this, id)
            m = this.m_size(1,1);
            E = eye(this.m_size);
            M = this.m_E * this.m_M;
            if (M(id,id) ~= 0)
                E(id+1:m, id) = -M(id + 1: m, id)/M(id,id);
            else
                E(id, id) = 0;
                E(id + 1, id + 1) = 0;
                E(id, id + 1) = 1;
                E(id + 1, id) = 1;
            end
            this.m_E = E * this.m_E;
            T = this.m_E * this.m_M;
            g = this.show(T);
        end
        function SetGraphPose(this, g)
            set(g, 'XData', this.m_P(1,:), 'YData', this.m_P(2, :));
        end
        function ShowMultiplication(this, A)
            T = A * this.m_M;
            this.show(T);
        end
    end
    methods (Static)
        function g = plot(G)
            figure
            axis off
            hold on;
            color = rand(1,3);
            g = plot(G, 'LineWidth', 1, 'NodeFontSize',15, 'EdgeColor', color);
        end
    end
end