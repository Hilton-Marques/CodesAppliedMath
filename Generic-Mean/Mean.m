classdef Mean < Scene
    properties
        m_pts;
        m_n = 0;
        m_fac = 0.002;
        m_handles = [];
        m_x = [1;1];
    end
    methods
        function this = Mean(pts)
            this@Scene("mean.gif")
            if nargin > 1
                this.m_pts = pts;
                this.m_n = size(pts, 2);
                this.Init();
            end
        end

        function GaussianSum(this)
            alpha = 0.85;
            m = 0;
            s = 0.5;
            n = 100;
            t = linspace(m - 3*s, m + 3*s, n);
            g = alpha * this.gaussian(m,s,t);
            this.showGaussian(g,t, this.m_red);
            m = 2;
            s = 0.5;
            t = linspace(m - 3*s, m + 3*s, n);
            g =  (1-alpha) * this.gaussian(m,s,t);
            this.showGaussian(g,t, this.m_blue);
            axis tight
        end

        function ShowDataGaussian(this, data, color)
            data = data(1,:); %Get just one dimensional projection
            m = mean(data);
            s = std(data);
            n = 100;
            t = linspace(m - 3*s, m + 3*s, n);
            g = this.gaussian(m,s,t);
            this.showGaussian(g,t, color);
        end

        function Solver(this)
            iter = 100000;
            %this.Show();
            this.get();
            for i = 1:iter
                g = -this.m_fac * this.Grad(this.m_x);
                this.m_x = this.m_x + g;
                if (mod(i,50) == 0)
                    %this.Show();
                    %this.get();
                end
            end
            %this.Show();
            %this.get();
            %this.save(repeat=true);
        end

        function Init(this)
            bb = this.setBB(this.m_pts);
            %axis(bb);
            %plot(this.m_pts(1,:), this.m_pts(2,:), 'o', 'MarkerFaceColor', this.m_red, 'MarkerSize', 8, 'MarkerEdgeColor','black');
        end

        function ShowPts(this, pts, color)
           plot(pts(1,:), pts(2,:), 'o', 'MarkerFaceColor', color, 'MarkerSize', 4, 'MarkerEdgeColor','black');
        end

        function Show(this)
            delete(this.m_handles);
            this.m_handles = [];
            
            for i = 1:this.m_n
                this.m_handles = [this.m_handles, ...
                                  this.Spring(this.m_pts(:,i), this.m_x, 0.01)];
            end
            this.m_handles = [this.m_handles, ...
                             plot(this.m_x(1,:), this.m_x(2,:), 'o', 'MarkerFaceColor', this.m_blue, 'MarkerSize', 10)];
        end

        function g = Grad(this, x)            
            g = [0;0];
            for i = 1:this.m_n
                g = g + (x - this.m_pts(:,i))/norm((x - this.m_pts(:,i)));
                %g = g + 2*(x - this.m_pts(:,i));
            end
            g = g/this.m_n;
        end

        function ShowCosts(this)
            fig = figure;
            I = 2.5;
            M = 1.5;
            n = 20;
            n2 = n - 2;
            %L2
            x1 = linspace(-I, I, n);
            u1 = x1.^2;            
            %L1 
            n = 5;
            x2 = linspace(-I,I,n);
            u2 = [abs(x2(1:n/2)), 0, flip(abs(x2(1:n/2)))];
            %Huber
            x3 = linspace(-M, M, n2);            
            u3 = x3.^2;
            r = linspace(M,I,5);
            r_u = M*(2*r - M);
            u3 = [r_u,u3,r_u];
            x3 = [-r,x3,r];
            %Corrupted gaussian
            n = 20;
            x5 = linspace(-I,I,n);
            a = 0.9;
            w = 3.5;
            u5 = -log(a * exp(-x5.^2) + (1 - a) * exp((-x5.^2)/(w^2))/w);
            
            %Ref
            n = 20;
            x4 = linspace(-M, M, n-2);            
            u4 = x4.^2;
            u4 = [M^2,u4,M^2];
            x4 = [-I,x4,I];

            this.plotErrors(x1,u1,'x','L(x)', line_width=1.0, color=this.m_blue, marker='o', marker_color=this.m_blue);
            this.plotErrors(x2,u2,'x','L(x)', line_width=1.0, color=this.m_red, marker='+', marker_color=this.m_red);
            this.plotErrors(x3,u3,'x','L(x)', line_width=1.0, color=this.m_green, marker='s', marker_color=this.m_green);
            this.plotErrors(x5,u5,'x','L(x)', line_width=1.0, color=this.m_marron, marker='*', marker_color='cyan');
            this.plotErrors(x4,u4,'x','L(x)');

            legends = ["L_2", "L_1", "\mathrm{Huber}", "\mathrm{mistura}","\mathrm{ideal}"];
            this.plotLegends(legends);
            axis equal;
            %close(fig);

        end
    end
    methods (Static)
        
    end
end