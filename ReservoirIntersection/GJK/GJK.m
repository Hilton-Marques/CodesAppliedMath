classdef GJK < handle
    properties
        m_convex_A
        m_convex_B
        m_tetra_A
        m_tetra_B
        m_tetra_M % Minkowiski tetrahedra
        m_tol = 9.9e-7
    end
    methods
        function this = GJK(convex_A,convex_B)
            this.m_convex_A = convex_A;
            this.m_convex_B = convex_B;
            this.m_tetra_A = Tetrahedra();
            this.m_tetra_B = Tetrahedra();
            this.m_tetra_M = Tetrahedra();
        end
        function [bool, p] = compute(this, planes)
            bool = true;
            p = [];

            d = (this.m_convex_A.centroide - this.m_convex_B.centroide);
            if dot(d,d) < this.m_tol^2
                bool = true;
                p = this.m_convex_A.centroide;
                return;
            end
            d = d/norm(d);
            %d0 = [0.4,0.7,0.6];
            vA = this.m_convex_A.supportfcn(d');
            vB = this.m_convex_B.supportfcn(-d');
            A = vA - vB;
            if dot(A,A) < this.m_tol^2
                p = vA;
                bool = true;
                return
            end
            this.m_tetra_M.append(A);
            this.m_tetra_A.append(vA);
            this.m_tetra_B.append(vB);
            d = -A;
            d = d/norm(d);
            iter = 0;
            while true
                d = d/norm(d);
                vA = this.m_convex_A.supportfcn(d');
                vB = this.m_convex_B.supportfcn(-d');
                P = vA - vB;
                a = dot(P,d)/norm(P)
                if dot(P,d) <= this.m_tol
                    bool = false;
                    return
                end
                this.m_tetra_M.append(P);
                this.m_tetra_A.append(vA);
                this.m_tetra_B.append(vB);
                [flag,d] = this.handleSimplex();
                if flag
                    bool = true;
                    break;
                end
                iter = iter + 1;
            end
            lam = this.m_tetra_M.getOriginBary();
            p = this.m_tetra_A.getBaryPoint(lam);
            values = zeros(12,1);
            for i = 1:size(planes,1)
                n = planes(i,1:3);
                d = planes(i,4);
                u = p + n * d;
                angle = dot(u,n);
                values(i) = angle;
                u = u / norm(u) ;
                angle = dot(u,n);
                if (angle >= 0.0)
                    bool = false;
                    %return
                end
            end
        end
        function [bool, d] = handleSimplex(this)
            if this.m_tetra_M.m_n == 2
                [bool,d] = this.m_tetra_M.GetPerpDir2Line();
                return
            elseif this.m_tetra_M.m_n == 3
                [bool,d] = this.m_tetra_M.GetPerpDir2Tri();
                return
            end
            [bool, d] = this.m_tetra_M.ContainOrigin(this.m_tetra_A,this.m_tetra_B);
        end
    end
end