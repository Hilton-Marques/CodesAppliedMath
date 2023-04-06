classdef Draw < handle
    properties
        m_f
        m_margin = 0.3
        m_points = [];
        m_format = ['%', sprintf('.%df', 2)]
    end
    methods
        function this = Draw(f,p)
            this.m_f = f;
            this.m_points(:,end+1) = p;
            hold on
            this.Init()
        end
        function Init(this)     
            this.plotContour(5);
            p1 = this.m_points(:,end);
            this.plotP(p1);
            text(p1(1),p1(2),strcat('[',num2str(p1(1),3),{' '},num2str(p1(2),3),']'),'HorizontalAlignment','right','VerticalAlignment','bottom');
            hold off
        end
        function Update(this,p2,d)
            format = this.m_format;
            fac2 = 1.0;
            p1 = this.m_points(:,end);
            this.m_points(:,end+1) = p2;
            hold on            
            this.plotP(p2);
            this.longline(p1,p2);
            L = norm(p2 - p1);
            v = (p2 - p1)/norm(p2 - p1);    
            if (L > 0.1)
            text(p2(1),p2(2),strcat('[',num2str(p2(1),format),{' '},num2str(p2(2),format),']'),'HorizontalAlignment','left','VerticalAlignment','top');
            end
            o = [-v(2),v(1)];
            if L < 1
                fac2 = L;
            end
            fac = 0.07;
            quiver(p1(1)+ fac*o(1),p1(2)+fac*o(2),fac2*v(1),fac2*v(2),'linewidth',2,'color','red');
            %text(v(1)+p1(1)+fac*o(1),v(2)+p1(2)+fac*o(2),strcat('[',num2str(d(1),3),{' '},num2str(d(2),3),']'),'HorizontalAlignment','left','VerticalAlignment','bottom');
            bb = this.GetBoundingBox();
            axis(bb);
        end
        function plotFinalP(this,p,color,iter)
            this.plotContour();
            format = this.m_format;
            this.plotP(p,color,8);
            if iter > 3
            text(p(1),p(2),strcat('[',num2str(p(1),format),{' '},num2str(p(2),format),']'),'HorizontalAlignment','right','VerticalAlignment','top');
            end
        end
        function Export(~,filename)
            exportgraphics(gca,strcat(filename,'.jpeg'),'Resolution',333);
        end
        function bb = GetBoundingBox(this, margin)
            if nargin == 1
                margin = this.m_margin;
            end
            x_min = min(this.m_points(1,:)) - margin;
            x_max = max(this.m_points(1,:)) + margin;
            y_min = min(this.m_points(2,:)) - margin;
            y_max = max(this.m_points(2,:)) + margin;            
            bb = [x_min,x_max,y_min,y_max];
        end
        function plotContour(this,margin)
            if nargin == 1
                margin = this.m_margin;
            end
            ax = this.GetBoundingBox(margin);
            n = 20;
            x = linspace(ax(1),ax(2),n);
            y = linspace(ax(3),ax(4),n);
            [X,Y] = meshgrid(x,y);
            Z = zeros(n,n);
            %Z = this.m_f(X,Y);
            for i = 1:n
                for j = 1:n
                    Z(i,j) = this.m_f([X(i,j),Y(i,j)]);
                end
            end
            contour(X,Y,Z,linspace(min(min(Z)),max(max(Z)),20))
        end
    end
    methods (Static)
        function plotP(x,color,size)
            if nargin == 1
                color = 'blue';
                size = 6;
            end
            if nargin == 2
                size = 6;
            end
            plot(x(1),x(2),'o','MarkerFaceColor',color,'MarkerSize',size);
        end
        function longline(p1,p2,width)
            if nargin ==2
                width = 2;
            end
            line([p1(1),p2(1)],[p1(2),p2(2)],'linewidth',width,'color','black');
        end
    end
end