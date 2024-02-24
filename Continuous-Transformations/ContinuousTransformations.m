%{
%Copyright (c) 2023 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Sunday, November 12th 2023, 11:48:51 pm
%Author: Hilton-Marques
%
%Description: A class to render continuous transformations in the plane
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

classdef ContinuousTransformations < Scene
    properties
        m_dt = 0.1;
        m_t = 1.0;
    end
    methods
        function this = ContinuousTransformations(filename)
            this = this@Scene(filename);
        end
        function IntegrateInfinitesimalTransformation(this, T)
            [X,Y] = this.CartesianCoordinates(10);
            this.setBB([[-5,5];[-5,5]]);
            n = 10; %iterations
            this.m_dt = this.m_t/n;
            m = 20;%fill in between
            marker_size_end = 5;
            marker_size = 1;
            lambda = exp(log(marker_size_end)/(n*m)) - 1;
            count = 0;
            plot(X,Y,'o', 'MarkerFaceColor', this.m_blue,'Color', this.m_blue, 'MarkerSize',marker_size);
            this.get();
            for i = 1:n
                V = T(X,Y);
                dX = V{1};
                dY = V{2};
                X_f = X + dX * this.m_dt;
                Y_f = Y + dY * this.m_dt;
                color_i = (this.m_blue + ((i-1)/n) * (this.m_red - this.m_blue));
                color_j = (this.m_blue + (i/n) * (this.m_red - this.m_blue));

                %fill in between
                for j = 1:m
                    X_j = X + (j/m)*(X_f - X);
                    Y_j = Y + (j/m)*(Y_f - Y);
                    color_k = color_i + (j/m) * (color_j - color_i);
                    marker_size = (1 + lambda) ^count;
                    plot(X_j,Y_j,'o', 'MarkerFaceColor', color_k,'Color', color_k,'MarkerSize',3);
                    count = count + 1;
                    this.get();
                end                
                X = X_f;
                Y = Y_f;
                marker_size = (1 + lambda) ^count;                
                plot(X,Y,'o', 'MarkerFaceColor', color_i,'Color', color_i, 'MarkerSize',3);                
                count = count + 1;
                this.get();                
            end
        end
        function LinearTransformation(this,T)
            x = linspace(-5,5,20);
            this.setBB([x;x]);
            [X,Y] = meshgrid(x);
            V = T(X,Y);
            X_f = V{1};
            Y_f = V{2};
            n = 100;
            marker_size_end = 5;
            marker_size = 1;
            lambda = exp(log(marker_size_end)/(n)) - 1;
            plot(X,Y,'o', 'MarkerFaceColor', this.m_blue,'Color', this.m_blue, 'MarkerSize',marker_size);
            %plot(X_f,Y_f,'o', 'MarkerFaceColor', this.m_red,'Color', this.m_red, 'MarkerSize',marker_size_end);
            %quiver(X,Y, X_f - X, Y_f - Y);
            this.get();
            for i = 1:n
                X_i = X + (i/n) * (X_f - X);
                Y_i = Y + (i/n) * (Y_f - Y);        
                marker_size = (1 + lambda) ^i;
                color_i = (this.m_blue + ((i-1)/n) * (this.m_red - this.m_blue));
                plot(X_i,Y_i,'o', 'MarkerFaceColor', color_i,'Color', color_i, 'MarkerSize',marker_size);
                %plot(X_i(end,1),Y_i(end,1),'o', 'MarkerFaceColor', 'green','Color', color_i, 'MarkerSize',marker_size);
                this.get(); 
            end
            %X_i(end,1)
            %Y_i(end,1)
        end
    end
    methods (Static)
        function [X,Y] = PolarCoordinates(n)
            r = linspace(-5,5,n);
            theta = linspace(0,2*pi,n);
            [R, Theta] = meshgrid(r,theta);
            X = R .* cos(Theta);
            Y = R .* sin(Theta);
        end
        function [X,Y] = CartesianCoordinates(n)
            x = linspace(-5,5,n);
            [X,Y] = meshgrid(x);
        end
    end
end

