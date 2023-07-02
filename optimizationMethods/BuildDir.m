classdef BuildDir < handle
    properties
        m_method_name
        m_f
        m_df
        m_ddf
        m_dim
        m_n_iter = 0
        m_d_powell    % Powell
        m_p           % Powell, BFGS
        m_p_end       % Powell 
        m_d           % FR, BFGS
        m_S           % BFGS
    end
    methods (Access = public)
        function this = BuildDir(method_name,f,p)
            this.m_method_name = method_name;
            this.m_f = f;
            this.m_dim = length(p);
            this.m_d_powell = eye(this.m_dim);
            this.m_p = p;
            this.m_S = eye(this.m_dim);
            %Init Derivatives
            if this.m_dim == 2
                this.m_df = this.Grad(f);
                this.m_ddf = this.Hessian(f);
                this.m_d = this.m_df(p);
            end
        end
        function out = GetDir(this,p)            
            switch this.m_method_name
                case 'Univariante'
                    out = this.GetDirUni();
                case 'Powell'
                    out = this.GetDirPow(p);
                case 'SD'
                    out = this.GetDirSD(p);
                case 'FR'
                    out = this.GetDirFR(p);
                case 'NR'
                    out = this.GetDirNR(p);
                case 'BFGS'
                    out = this.GetDirBFGS(p);
                case 'dogleg'
                    out = this.GetDirDogLeg(p);
            end
            this.m_n_iter = this.m_n_iter + 1;
        end
    end
    methods (Access = private)
        function d = GetDirUni(this)
            d = zeros(this.m_dim,1);
            pose = mod(this.m_n_iter,this.m_dim) + 1;
            d(pose,1) = 1.0;
        end
        function d = GetDirPow(this,p)
            pose = mod(this.m_n_iter,this.m_dim + 1) + 1;
            if pose == 1
              this.m_p = p;
            elseif pose == this.m_dim + 1
              this.m_p_end = p;
            end            
            if (pose == this.m_dim + 1)
                d = this.m_p_end - this.m_p;
                this.m_d_powell(:,1) = [];
                this.m_d_powell(:,this.m_dim) = d;
                %this.m_p = this.m_p_end;
            else
                d = this.m_d_powell(:,pose);
            end
            % check end of a cycle
            if mod(this.m_n_iter + 1,(this.m_dim + 1)^2) == 0
                this.m_d_powell = eye(this.m_dim);
                %this.m_p = this.m_p_end;
            end
        end
        function d = GetDirSD(this,p)
            d = -this.m_df(p);
        end
        function d = GetDirFR(this,p)
            d = this.m_d;
            if (this.m_n_iter == 0)
                d = -d;
                this.m_d = d;
                return;
            end
            grad = this.m_df(p);
            beta = dot(grad,grad)/dot(d,d);
            d = -grad + beta*d;
            this.m_d = d;
        end
        function d = GetDirNR(this,p)
            ddf = this.m_ddf(p);
            df = this.m_df(p);
            d = -ddf\df;
        end
        function d = GetDirBFGS(this,p)
            d = -this.m_d;
            if (this.m_n_iter == 0)
                return;
            end
            df = this.m_df(p);
            s = p - this.m_p;
            y = df - this.m_d;
            Sk = this.m_S;
            A = (s'*y + y'*Sk*y)*(s*s')/(s'*y)^2;
            B = (Sk*y*s' + s*y'*Sk)/(s'*y);
            S = Sk + A - B;
            d = -S*df;
            this.m_p = p;
            this.m_d = df;
            this.m_S = S;            
        end
        function d = GetDirDogLeg(p)
        end
    end
    methods (Static)
        function out = Grad(f)
            h = 1e-5;
            %out = @(x,y) [imag(f(x+ 1i*h,y))/h ; imag(f(x,y + 1i*h))/h];
            out = @(x) [imag(f([x(1)+ 1i*h,x(2)]))/h ; imag(f([x(1),x(2) + 1i*h]))/h];
        end
        function out = Hessian(f)
            h = 1e-5;
            %             out = @(x,y) [[imag(f(x + 1i*h + h,y) - f(x + 1i*h - h,y))/(2*h^2),imag(f(x + 1i*h,y + h) - f(x + 1i*h ,y-h))/(2*h^2)];...
            %                [imag(f(x+h,y + 1i*h) - f(x-h,y + 1i*h))/(2*h^2) ,imag(f(x,y + 1i*h + h) - f(x,y + 1i*h - h))/(2*h^2)]];
            out = @(x) [[imag(f([x(1) + 1i*h + h,x(2)]) - f([x(1) + 1i*h - h,x(2)]))/(2*h^2), ... 
                         imag(f([x(1) + 1i*h,x(2) + h]) - f([x(1) + 1i*h ,x(2) - h]))/(2*h^2)];...
                         [imag(f([x(1) + h,x(2) + 1i*h]) - f([x(1) - h,x(2) + 1i*h]))/(2*h^2) ,...
                         imag(f([x(1),x(2) + 1i*h + h]) - f([x(1),x(2) + 1i*h - h]))/(2*h^2)]];
        end
    end
    
end