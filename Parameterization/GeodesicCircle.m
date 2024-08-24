classdef GeodesicCircle < Scene
	properties
		
	end

	methods
		function this = GeodesicCircle()
			this = this@Scene("GeodesicCircle", "S2");
		end

		function Solver(this,P,Q)
            view(65,30)
            axis tight
			s = S2();
            v = s.logMap(P,Q);
            t = linspace(0,1,20);
            g = s.geodesic(v,P,t);
            n = size(g,1);
            
            fig = figure;
            ax = axes(fig);
            set(gcf,'color','white');
            xlabel(ax,"$\theta$","Interpreter","latex");            
            ylabel(ax,"$\phi$","Interpreter","latex");
            axis equal

            hold on
            axis([-pi,pi,-pi/2,pi/2]);
            plot3(this.m_axes,Q(1),Q(2),Q(3),'o', 'MarkerSize',5, 'MarkerFaceColor',this.m_green);
            color = this.m_red;
            for i = 1:n
                pts = s.GeodesicCircle(g(i,:), 0.25);
                this.drawPolygon(pts',this.m_axes,color);
                %plot3(this.m_axes, pts(:,1), pts(:,2), pts(:,3),'o','MarkerFaceColor','red','MarkerSize',10);                
                [theta, phi] = this.GetLongLat(pts);           
                
                this.drawPolygon([theta,phi]',ax,color);                
                color = (this.m_red + (i/n)*(this.m_blue - this.m_red));
                %plot(ax, theta, phi, '-','color',color);
            end
            exportgraphics(ax, "r2_projection.jpeg","Resolution",300);
            %this.exportFrame();
    		
		end
    end

    methods (Static)
        function [theta, phi] = GetLongLat(pts)
			[theta, phi] = cart2sph(pts(:,1), pts(:, 2), pts(:,3));
        end
    end

end