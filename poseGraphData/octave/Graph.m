classdef Graph < handle
    properties
        m_poses_0
        m_poses
        m_edges
        m_drawer
        m_G
        m_n
        m_m
    end
    methods
        function this = Graph(poses,edges)
            n = size(poses,1);
            m = size(edges,1);
            this.m_poses = poses';
            this.m_poses_0 = poses';
            this.m_edges = edges;
            this.m_drawer = Drawer();
            this.m_n = n;
            this.m_m = m;
            %this.ErrorDeformation()
            %this.Show();
            %this.exportFrame('first');
            %pause(1);
        end
        function Show(this)
            this.m_drawer.Reset();
            len = -1;
            for i = 1:this.m_m
                pi = this.m_poses(:,this.m_edges(i,1));
                pj = this.m_poses(:,this.m_edges(i,2));
                li = dot(pj - pi,pj - pi);
                if (li > len)
                    len = li;
                end
                e = pj - pi;
                angle = atan2(e(2),e(1));
                %this.m_drawer.DrawEdges(pi,pj);
                this.m_drawer.DrawSprings(pi,pj);
                if (i == 1)
                    this.m_drawer.DrawNode([pi;angle],this.m_drawer.m_red);
                else
                    this.m_drawer.DrawNode([pi;angle],this.m_drawer.m_blue);
                end
                %this.m_drawer.update();
                %pause(0.5);
            end
            this.m_drawer.update();
        end
        function UpdatePose(this,pose)
            this.m_poses = pose';
        end
        function star = GetStar(this,i)
            ids = 1:this.m_m;
            star = this.m_edges(i == this.m_edges(:,1),2);
        end
        function Save(this)
            this.m_drawer.save('show_Map_tot.gif');
        end
        function ErrorDeformation(this)
            axis([-0.5,2.2,-0.5,1.5]);
            c = [0;0];
            p0 = [1.0;0.0];
            [c,d] = this.trajectory(c,p0);
            n = length(c);
            p0(1) = p0(1)  - 0.1;
            %this.m_drawer.origin(p0,d(:,1));
            plot(-0.1,1.0,'x','markerfacecolor','green','markersize',15,'linewidth',2,'color','green');

            h = this.m_drawer.DrawNode([c(:,1);atan2(d(2,1),d(1,1))],this.m_drawer.m_blue);
            this.m_drawer.update();
            for i = 2:18
                h = this.m_drawer.DrawNode([c(:,i);atan2(d(2,i),d(1,i))],this.m_drawer.m_blue);
                plot(c(1,i),c(2,i),'o','color','black','markerfacecolor','black','markersize',1);
                this.m_drawer.update();
                delete(h);
            end
            this.m_drawer.DrawNode([c(:,i);atan2(d(2,i),d(1,i))],this.m_drawer.m_red);
            this.m_drawer.update();
            %this.m_drawer.save("error.gif");
            %plot(0.2,0.8,'x','markerfacecolor','black','markersize',15,'linewidth',2,'color','black');
            this.m_drawer.save("error_2.gif");
            this.drawArc(c(:,3),c(:,i),-1,"$T^{-}_{x_{j}}T_{x_{k}}$",'top');

            keyboard;
            this.drawArc([0;0],c(:,1),1,"$T_{x_i}$",'top');
            this.m_drawer.export("Txi");
            this.drawArc([0;0],c(:,i),-1,"$T_{x_{i+1}}$",'top');
            this.m_drawer.export("Txi2");
            this.drawArc(c(:,3),c(:,i),1,"$u_{x_{i};x_{i+1}}$",'bottom');
            this.m_drawer.export("odometry");
        end
        function exportFrame(this,filename)
            this.m_drawer.export(filename)
        end
        
    end
    methods(Static)
        function [c,d] = trajectory(c,p0)
            u = p0 - c;
            L = norm(u);
            u = u/L;
            v = [-u(2);u(1)];
            u = L * 0.9 * u;
            v = (L * 0.4 * v);
            theta = linspace(0,2*pi,100);
            c = cos(theta) .* u  + sin(theta) .* v +c;
            d = diff(c')';
            d(:,end+1) = d(:,1);
            d = d ./ vecnorm(d);
        end
        function c = drawArc(p0,p1,sign,label,pose)
            u = p0 - p1;
            L = norm(u);
            u = u/L;
            v = sign*[-u(2);u(1)];            
            theta = linspace(0,pi - pi/12,100);
            u = L * 0.5 * u;
            v = (L * 0.25 * v);
            c = cos(theta) .* u  + sin(theta) .* v + p0 - u;
            plot(c(1,:),c(2,:),'linewidth',1,'color','black');
            pf = c(:,end);
            d = c(:,end) - c(:,end-1);
            angle = atan2(d(2),d(1));
            delta = pi/6;
            r = L/8;
            left = pf - r*[cos(angle - delta);sin(angle - delta)];
            right = pf - r*[cos(angle + delta);sin(angle + delta)];
            line([pf(1),left(1)],[pf(2),left(2)],'Color','black','linewidth',1);
            line([pf(1),right(1)],[pf(2),right(2)],'Color','black','linewidth',1);
            text(c(1,50),c(2,50),label,'interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment',pose);
            
        end
    end
end