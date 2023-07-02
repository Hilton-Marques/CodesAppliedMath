classdef FindMinLinSearch < handle
    properties
        m_dim
        m_tol = 10^-5
        m_alpha = 0.001;
        m_max_iter = 500;
        m_lin_searcher
        m_dir_builder
        m_drawer
        m_p0
        m_last_iter = -1
    end
    methods
        function [this,p_min] = FindMinLinSearch(f,p,lin_method,dir_method)
            this.m_p0 = p;
            this.m_lin_searcher = LinearSearch(this.m_tol,this.m_alpha,lin_method,f);
            this.m_dir_builder  = BuildDir(dir_method,f,p);
            this.m_dim = length(p);
            if this.m_dim == 2
                this.m_drawer = Draw(f,p);
            end
            p_min = this.Solve(p);
        end
        function p_min = Solve(this,p)
            iter = 0;
            tic
            while iter < this.m_max_iter
                d = this.m_dir_builder.GetDir(p);
                new_p = this.m_lin_searcher.GetMinPoint(p,d);
                this.UpdateScene(new_p,d);
                if this.IsFinished(new_p,p)
                    iter = iter + 1;
                    break;
                end
                p = new_p;
                iter = iter + 1;
            end
            toc
            this.m_last_iter = iter
            p_min = new_p;
            this.UpdateScene(new_p,d);
            %% Export image
            %this.m_drawer.Export(strcat(this.m_dir_builder.m_method_name,num2str((this.m_p0')),'iter=',num2str(iter)))
        end
        function bool = IsFinished(this,new_p,p)
            bool = false;
            x_diff = norm(new_p - p)/norm(new_p);
            if this.m_dim == 2
                grad_p = norm(this.m_dir_builder.m_df(new_p));
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
end