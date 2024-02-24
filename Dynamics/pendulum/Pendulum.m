%{
%Copyright (c) 2023 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Saturday, September 23rd 2023, 12:37:37 pm
%Author: Hilton-Marques
%
%Description: A simple clase to study the pendulum dynamics
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

classdef Pendulum < Scene
    properties
        m_L = 1.0;
        m_M = 1.0;
        m_theta = 0.0;
        m_v = 0.0;
        m_h = [];
        m_bb = [-1,1,-1.2,0.2];
    end
    methods (Access=public)
        function this = Pendulum(length, mass)
            this = this@Scene("pendulum.gif");
            if nargin > 0
                this.m_L = length;
                this.m_M = mass;
            end
            this.Show();
            axis(this.m_bb);
        end

        function Trajectory(this, theta_fcn)
            n = 50;
            t = linspace(0, 1, n);
            thetas = theta_fcn(t);
            for i = 1:n
                ti = i/n;
                this.m_theta = thetas(i);
                this.m_v = this.Dtheta(theta_fcn, ti);
                this.Show();
                traj = this.Coord(theta_fcn(linspace(0, ti, i)));                
                h = this.drawPolygon(traj);
                this.m_h = [this.m_h, h];
                if (mod(i, 5) == 0 || i == 1)
                    this.DrawCircle();
                    this.DrawVelocity();
                end
                this.get();
                %pause(0.1);
            end
        end

        function Show(this)
            delete(this.m_h);            
            c = this.Coord();            
            %h1 = line([0,c(1)], [0,c(2)], 'linewidth', 2, 'color', 'black');
            h2 = this.DrawCircle();
            h3 = this.DrawVelocity();
            this.m_h = [h2, h3];
        end
    end
    methods (Access=private)
        function p = Coord(this, theta)
            if (nargin == 1)
                theta = this.m_theta;
            end
            y = -this.m_L * cos(theta);
            x = this.m_L * sin(theta);
            p = [x;y];
        end

        function h = DrawCircle(this)
            h = this.drawCircle(this.Coord(), 0.05, [0.5, 0.5, 0.5]);
        end
        
        function h = DrawVelocity(this)
            c = this.Coord();
            t = 0.07 * this.m_v * [-c(2); c(1)];
            h = quiver(c(1), c(2), t(1), t(2), 'color', this.m_red, 'linewidth', 2.0);
        end

        function v = Dtheta(this, theta_fcn, t)
            h = 1e-6;
            dtheta = (theta_fcn(t + h) - theta_fcn(t))/h;
            v = dtheta * this.m_L;
        end
    end
    methods (Static)
        
    end
end