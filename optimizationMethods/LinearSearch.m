classdef LinearSearch < handle
    properties (Constant)
        m_all_methods = ["Bisseção", "Seção Áurea"];
    end
    properties
        m_tol
        m_d
        m_alpha
        m_lower
        m_upper
        m_method_name % método de buscar Bisseção, Áurea
        m_f % function
        
    end
    methods
        function this = LinearSearch(tol,alpha,searcher_method, f)
            this.m_tol = tol;
            this.m_alpha = alpha;
            this.m_method_name = searcher_method;
            this.m_f = f;
        end
        function out = GetMinPoint(this,p,d)
            this.PassoConstante(p,d);
            this.m_d = this.m_d/norm(this.m_d);
            switch this.m_method_name
                case 'Bisseção'
                    out = this.GetPointBissec();
                case 'Seção Áurea'
                    out = this.SecAurea();
            end
        end        
    end
    methods (Static)
        function bool = IsValidMethod(method)
            bool = false;
            for name = LinearSearch.m_all_methods
                if method == name
                    bool = true;
                    return;
                end
            end
        end
    end
    methods (Access = private)
        function xi = PassoConstante(this, p, d)
            this.m_d = d;
            f0 = this.m_f(p);
            xi = p + this.m_tol*d;
            if (this.m_f(xi) >= f0)
                this.m_d = -this.m_d;
            end
            xi = p + this.m_alpha*this.m_d;
            while (this.m_f(xi) < f0)
                f0 = this.m_f(xi);
                xi = xi + this.m_alpha*this.m_d;
            end
            this.m_lower = xi - this.m_alpha*this.m_d;
            this.m_upper = this.m_lower + this.m_alpha*this.m_d;
        end
        function out = GetPointBissec(this)
            p1 = this.m_lower;
            p2 = this.m_upper;
            D = norm(p2 - p1);
            D = D/2;
            M = p1 + D*this.m_d;
            while (D > this.m_tol)                
                M_l = M - this.m_tol*this.m_d;
                M_u = M + this.m_tol*this.m_d;
                L = this.m_f(M_l);
                U = this.m_f(M_u);
                if ( U < L)
                    p1 = M;
                end
                D = D/2;
                M = p1 + D*this.m_d;
            end
            out = M;
        end
        function out = SecAurea(this)
            p1 = this.m_lower;
            p2 = this.m_upper;
            d = this.m_d;
            Ra = (sqrt(5) - 1)*0.5;
            D = norm(p2 - p1);
            M_u = p1 + Ra*D*d;
            M_l = p1 + (1 - Ra)*D*d;
            L = this.m_f(M_l);
            U = this.m_f(M_u);
            while (D > this.m_tol)
                D = Ra*D;
                if ( U < L)
                    M_l = M_u;
                    p1 = M_l;
                    L = U;
                    M_u = p1 + Ra*D*d;
                    U = this.m_f(M_u);
                else
                    M_u = M_l;
                    U = L;
                    M_l = p1 + (1 - Ra)*D*d;
                    L = this.m_f(M_l);
                end
            end
            out = (M_u+M_l)/2;
        end
    end
end