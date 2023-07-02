classdef Elipse < handle 
    properties 
        m_a;
        m_b;
        m_center;
        m_l
        m_h
    end
    methods
        function this = Elipse(center,a,b)
            this.m_center = center;
            this.m_a = a;
            this.m_b = b;   
            r = (a + b)*0.5;
            this.m_l = r;
        end
        function show(this)
            teta = linspace(0,2*pi,30);
            x = this.m_center(1) + cos(teta)*this.m_a;
            y = this.m_center(2) + sin(teta)*this.m_b;
            plot(x,y,'linewidth',1,'color','black');
            fill(x,y,'yellow','FaceAlpha',0.5);
            [x_text,y_text] = this.getTextPose(pi/4);            
            this.m_h = text(x_text,y_text,'$B,A$','FontSize',10,'FontWeight','bold','Interpreter','latex');
            this.Partition();
        end
        function [x,y] = getTextPose(this,teta)
            p = this.m_center + this.m_l * [cos(teta),sin(teta)];
            x = p(1);
            y = p(2);
        end
        function Partition(this)
            h = 2*pi/4;
            poses_l = ["right","left","left","right"];
            poses_t = ["top","top","bottom","bottom"];

            delete(this.m_h)
            for i = 1:4
                teta = linspace((i-1)*h,i*h,30);
                x = this.m_center(1) + cos(teta)*this.m_a;
                y = this.m_center(2) + sin(teta)*this.m_b;
                x(end+1) = this.m_center(1);
                y(end+1) = this.m_center(2);
                fill(x,y,rand(1,3),'FaceAlpha',0.5);
                str = strcat("$B|A_{",num2str(i),"}$");
                [x,y] = this.getTextPose((i-1)*pi/2 + pi/4);
                text(x,y,str,'FontSize',10,'FontWeight','bold','Interpreter','latex','VerticalAlignment',poses_t(i),'HorizontalAlignment',poses_l(i));
            end

        end
    end
end