classdef Drawer < handle
    properties
        m_h = [];
        m_e = [];
        m_figure
        im = {};
        m_red = [0.9176 0.2627 0.2078];
        m_blue = [0.2588 0.5216 0.9569];
        filename = 'graph_force'
        m_r = 0.1
    end
    methods
        function this = Drawer()
            this.m_figure = figure;
            hold on
            axis equal
            axis off
            set(gcf,'color','w');
        end
        function update(this)
            frame = getframe(this.m_figure);
            this.im{end+1} = frame2im(frame);
        end
        function origin(this,p0,dir)
            o = [0,0];
            d = [0,1];
            if nargin > 1
                o = p0
                d = dir
            end
           d_ortho = [d(2),-d(1)];
            plot(o(1),o(2),'x','markersize',5);
            quiver(o(1),o(2),d(1),d(2),'color','black');
            quiver(o(1),o(2),d_ortho(1),d_ortho(2),'color','black');
        end
        function SetRadius(this,r)
            this.m_r = r;
        end
        function h = DrawRobot(this,x)
            r = 0.25;
            x = pose(1:2,:);
            d = r*[cos(pose(3,:));sin(pose(3,:))];
            h = plot(this.m_axis(ax),x(1,:),x(2,:),'o','MarkerFaceColor',color,'markersize',5,'markeredgecolor','black');
            h(end+1) = quiver(this.m_axis(ax),x(1,:),x(2,:),d(1,:),d(2,:),'linewidth',1,'color','black','autoscale','off','ShowArrowHead','off');
            teta = linspace(0,2*pi,20);
            n = size(pose,2);
            if n == 2
                color = [[0,0,1];0.85*[0,1,0]];
            end
            for i = 1:n
                x = pose(1:2,i);
                c = x + r*[cos(teta);sin(teta)];
                h(end+1) = plot(this.m_axis(ax),c(1,:),c(2,:),'linewidth',2,'color',color(i,:));
            end
        end
        function h = DrawNode(this, pose,color)
            r = 0.1;
            x = pose(1:2);
            theta = linspace(0,2*pi,30);
            circle = x + r*[cos(theta);sin(theta)];
            %this.m_h(end+1) = plot(x(1,:),x(2,:),'Color','red');
            h = fill(circle(1,:),circle(2,:),color);
            this.m_h(end+1) = h;            
            d = r*[cos(pose(3,:));sin(pose(3,:))];
            h2 = quiver(x(1,:),x(2,:),d(1,:),d(2,:),'linewidth',1,'color','black','autoscale','off','ShowArrowHead','off');
            return;
            d_ortho = [-d(2);d(1)];
            p0 = x + 5*d - 2.5*d_ortho;
            p1 = x + 5*d + 2.5*d_ortho;
            h_3 = line([p0(1),p1(1)],[p0(2),p1(2)],'linewidth',2,'color','black');
            h_4 = line([p0(1),x(1)],[p0(2),x(2)],'linewidth',0.5,'color','black','LineStyle','--');
            h_5 = line([p1(1),x(1)],[p1(2),x(2)],'linewidth',0.5,'color','black','LineStyle','--');
            h(end+1) = h2;
            h(end+1) = h_3;
            h(end+1) = h_4;
            h(end+1) = h_5;
            h(end+1) = h2;            
        end
        function DrawEdges(this,p1,p2)
            this.m_e(end+1) = line([p1(1),p2(1)],[p1(2),p2(2)],'Color','black');
        end
        function export(this,filename)
            exportgraphics(this.m_figure,strcat(filename,'.jpeg'),'Resolution',300);
        end
        function DrawSprings(this,p0,p1)
            n = 20;
            l = norm(p1 - p0);
            lam = 0.05/l;
            u = p1 - p0;
            ortho = [-u(2);u(1)];

            h = l/n;
            pi = p0;
            u = u/l;
            for i = 1:n-1
                pj = p0 + i*h*u + lam*ortho*(-1)^i;
                %handles(end+1) = line(this.m_axis(this.m_ax),[pi(1),pj(1)],[pi(2),pj(2)]);  
                this.m_e(end+1) = line([pi(1),pj(1)],[pi(2),pj(2)],'Color','black');                
                pi = pj;
            end
            pj = p1;
            this.m_e(end+1) = line([pi(1),pj(1)],[pi(2),pj(2)],'Color','black');
            
           % this.m_e(end+1) = line([p1(1),p2(1)],[p1(2),p2(2)],'Color','black');
        end
        function save(this,filename)
            for idx = 1:length(this.im)
                [A,map] = rgb2ind(this.im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
                elseif (idx > 1) && (idx < length(this.im))
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.08);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',5);
                end
            end
        end
        function Reset(this)
            delete(this.m_h);
            delete(this.m_e);
            this.m_h = [];
            this.m_e = [];
        end
    end
end