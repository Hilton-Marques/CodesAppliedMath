classdef FindMinTrustRegion < handle
    properties (Constant)
        m_tol = 1e-4;
        m_min_rel = 0.15;
        m_max_radius = 100;
        m_max_iter = 100;
    end
    properties (Access = public)
        m_method_name;
        m_f;
        m_g;
        m_h;
        m_aprox_q_solver
        m_radius = 1.0;
        m_drawer;
        m_dim;
        m_last_iter = -1
    end
    methods (Access = public)
        function [this, p_min] = FindMinTrustRegion(f, p,method_name)
            this.m_method_name = method_name;
            this.m_f = f;
            this.m_g = this.Grad(f);
            this.m_h = this.Hessian(f);
            this.m_dim = length(p);
            if this.m_dim == 2
                this.m_drawer = Draw(f,p);
            end
            this.SetMethod();
            p_min = this.Solve(p);
        end
        function p_min = Solve(this, p)
            iter = 0;
            while (iter < this.m_max_iter)
                u = this.m_g(p);
                H = this.m_h(p);
                v = - H \u;
                quad = dot(u,H*u);
                d = this.m_aprox_q_solver(u, v, quad);
                rel = (this.m_f(p + d) - this.m_f(p)) /(dot(u,d) + 0.5*(dot(d,H*d)));
                norm_d = norm(d);
                while (rel < this.m_min_rel)
                    if rel < 0.25
                        this.m_radius = 0.25 * norm_d;
                    end
                    d = this.m_aprox_q_solver(u, v, quad);
                    rel = (this.m_f(p + d) - this.m_f(p)) /(dot(u,d) + 0.5*(dot(d,H*d)));
                    norm_d = norm(d);
                end
                if rel > 0.75 && norm_d == this.m_radius %abs(norm_d - this.m_radius) < 1e-8
                    this.m_radius = min (2*this.m_radius, this.m_max_radius);
                end
                new_p = p + d;
                this.UpdateScene(new_p,d);
                if this.IsFinished(new_p,p)
                    iter = iter + 1;
                    break;
                end
                p = new_p;
                iter = iter + 1;
            end
            this.m_last_iter = iter
            p_min = new_p;
            this.UpdateScene(new_p,v);
        end
        function SetMethod(this)
            switch this.m_method_name
                case 'dogleg'
                    this.m_aprox_q_solver = @this.DogLeg;
                case 'Levenberg'
                    this.m_aprox_q_solver = @this.Levenberg;
            end 
        end
    end
    methods (Access = private)
        function d = DogLeg(this, b, d, quad)           
            dot_r = dot(d,d);
            radius_sqr = this.m_radius * this.m_radius;
            if dot_r <= radius_sqr
                %means that the trust region is big enough and encompass
                %the minimum
                return;
            end
            u = (-dot(b,b) / quad) * b; % cauchy point minimun distance 
                                                % along the line b
            dot_u = dot(u,u);
            if (dot_u >= radius_sqr)
                %means that the trust region is very small and the steppest
                %descent is enough
                d = (this.m_radius / sqrt(dot_u)) * u;
                return;
            end
            v = d - u;
            c = dot(v,v);
            a = dot_u;
            b = dot(v,u);
            delta = b^2 - a *c + c*radius_sqr;  % (a - 2b + c -d) + 2(b-c)t + ct^2
            sqr_delta = sqrt(delta);
            t1 = ((c - b) + sqr_delta) / c;
            t2 = ((c - b) - sqr_delta) / c;
            if t1 >= 1 && t1 <=2 
                t = t1;
            elseif t2 >= 1 && t2 <=2 
                t = t2;
            else
                disp('error');
            end
            d = u + (t - 1) * v;
        end
        function bool = IsFinished(this,new_p,p)
            bool = false;
            x_diff = norm(new_p - p)/norm(new_p);
            if this.m_dim == 2
                grad_p = norm(this.m_g(new_p));
                if (grad_p < this.m_tol*10 || x_diff < this.m_tol)
                    bool = true;
                end
            end
        end
        function UpdateScene(this,new_p,d)
            if this.m_dim == 2
                this.m_drawer.Update(new_p,d);
                if this.m_last_iter ~= -1
                    this.m_drawer.plotFinalP(new_p,'cyan',this.m_last_iter);
                end
            end
        end
    end
    methods (Static)
        function out = Grad(f)
            h = 1e-6;
            %out = @(x,y) [imag(f(x+ 1i*h,y))/h ; imag(f(x,y + 1i*h))/h];
            out = @(x) [imag(f([x(1)+ 1i*h,x(2)]))/h ; imag(f([x(1),x(2) + 1i*h]))/h];
        end
        function out = Hessian(f)
            h = 1e-6;
            %             out = @(x,y) [[imag(f(x + 1i*h + h,y) - f(x + 1i*h - h,y))/(2*h^2),imag(f(x + 1i*h,y + h) - f(x + 1i*h ,y-h))/(2*h^2)];...
            %                [imag(f(x+h,y + 1i*h) - f(x-h,y + 1i*h))/(2*h^2) ,imag(f(x,y + 1i*h + h) - f(x,y + 1i*h - h))/(2*h^2)]];
            out = @(x) [[imag(f([x(1) + 1i*h + h,x(2)]) - f([x(1) + 1i*h - h,x(2)]))/(2*h^2), ...
                imag(f([x(1) + 1i*h,x(2) + h]) - f([x(1) + 1i*h ,x(2) - h]))/(2*h^2)];...
                [imag(f([x(1) + h,x(2) + 1i*h]) - f([x(1) - h,x(2) + 1i*h]))/(2*h^2) ,...
                imag(f([x(1),x(2) + 1i*h + h]) - f([x(1),x(2) + 1i*h - h]))/(2*h^2)]];
        end
    end
end