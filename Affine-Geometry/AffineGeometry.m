classdef AffineGeometry < Scene
    properties
    end
    methods
        function this = AffineGeometry(filename)
            this = this@Scene(filename);
        end
        function TransformTriangle(this, p, q, r)
            this.setBB([p,q,r,[0;0]]);
            A = [q - p, r - p];
            n = 10;
            this.ShowTriangle([0;0],[1;0],[0;1]);    
            this.get();
            for i = 1:n
                A_i = this.Interp3DMatrix(A, i/n);
                %A_i = eye(2) + (i/n)*(A - eye(2));               
                p_i = (i/n)*p;
                q_i = p_i + A_i*[1;0];
                r_i = p_i + A_i*[0;1];
                color_i = (this.m_blue + (i/n) * (this.m_red - this.m_blue));
                h = this.ShowTriangle(p_i, q_i, r_i, color_i);
                norm(cross([q_i - p_i;0], [(r_i - p_i);0]));
                this.get();
            end
        end

        function TransformPolygon(this,A,O)
            if (size(A,2) == 2)
                A(3,3) = 1.0;
                O(3,:) = 0.0;
            end
            n = 20;
            h = this.FillPolygon(O, this.m_blue);
            this.get();
            for i = 1:n
                A_i = this.LinearInterp3DMatrix(A, i/n);                
                O_i = A_i * O;
                color_i = (this.m_blue + (i/n) * (this.m_red - this.m_blue));
                %delete(h);
                h = this.FillPolygon(O_i,color_i);
                this.get();
            end
        end
    end
    methods (Static)

    end
    
end