classdef Beam < handle
    properties
        m_h
        m_l
        m_k = 0;
        m_theta = 0;
        m_side
        m_pts
        m_margin = 2;
        m_handles = []
        m_blue = [0.2588 0.5216 0.9569];
    end
    methods
        function this = Beam(h,l)
            this.m_h = h;
            this.m_l = l;
            this.m_pts = [[0,0];[this.m_l,0.0];[this.m_l, this.m_h];[0,this.m_h]]';
            this.m_side = [0;this.m_h];
        end
        function updateCurvature(this,theta)
            %theta*R = L, k = theta/L
            R = this.m_l /(2*sin(theta));
            this.m_k = 1/R;
            this.m_theta = theta;
        end
        function bb = getBB(this)
            med = mean(this.m_pts');
            x_min = med(1) - (this.m_h/2) - this.m_margin;
            x_max = med(1) + (this.m_h/2) + this.m_margin;
            y_max = med(2) + (this.m_h/2) + this.m_margin/2;
            y_min = med(2) - (this.m_h/2) - this.m_margin;
            bb = [x_min,x_max,y_min,y_max];
        end
        function plot(this)
            region = [];
            if (this.m_k == 0)
                for i = 1:4
                    l = this.plotLine(this.m_pts(:,i),this.m_pts(:,mod(i,4)+1));  
                    region = [region,l];
                end
            else
                rot = this.Rot(this.m_theta);
                left = rot*this.m_side;
                right = rot'*this.m_side;
                this.m_pts(:,4) = left;
                this.m_pts(:,3) = this.m_pts(:,2) + right;
                radius = (1/this.m_k);
                center = -(left/this.m_h)*radius;
                %check = (this.m_pts(:,2) - (right/this.m_h)*radius) - center;
                theta1 = atan2(left(2),left(1));
                theta2 = theta1 - 2*this.m_theta;
                c1 = this.drawArc(center,radius,theta1,theta2);
                region = [region,c1];
                c2 = this.plotLine(this.m_pts(:,2),this.m_pts(:,3));
                region = [region,c2];
                c3 = this.drawArc(center,radius+this.m_h, theta2,theta1);
                region = [region,c3];
                c4 = this.plotLine(this.m_pts(:,4),this.m_pts(:,1));
                region = [region,c4];
                h1 = line([center(1),this.m_pts(1,1)],[center(2),this.m_pts(2,1)],'linewidth',1,'linestyle','--','color','black');
                h2 = line([center(1),this.m_pts(1,2)],[center(2),this.m_pts(2,2)],'linewidth',1,'linestyle','--','color','black');
                h3 = plot(center(1),center(2),'o','markersize',5,'markerfacecolor','black');
                h4 = text(center(1),center(2),'C','interpreter','latex','verticalalignment','top','horizontalalignment','left');
                this.m_handles(end+1) = h1;
                this.m_handles(end+1) = h2;
                this.m_handles(end+1) = h3;
                this.m_handles(end+1) = h4;
                %disp(this.length(c1));
            end
            h = fill(region(1,:),region(2,:),this.m_blue);
            this.m_handles(end+1) = h;
        end
        function clean(this)
            delete(this.m_handles);
            this.m_handles = [];
        end
        function arc = drawArc(this,center, radius, theta1,theta2)
            n = 100;
            theta = linspace(theta1,theta2,n);
            arc = center + radius*[cos(theta);sin(theta)];
            h = plot(arc(1,:),arc(2,:),'linewidth',1,'color','blue');
            this.m_handles(end+1) = h;
        end
        function l = plotLine(this, p0,p1)
            l = [p0,p1];
            h = line([p0(1,1),p1(1,1)],[p0(2,1),p1(2,1)],'linewidth',1,'color','blue');
            this.m_handles(end+1) = h;
        end
    end
    methods (Static)
        function rot = Rot(theta)
            rot = [[cos(theta),-sin(theta)];[sin(theta),cos(theta)]];
        end
        function s = length(curve)
            n = size(curve,2);
            s = 0;
            for i = 1:n-1
                pi = curve(:,i);
                pj = curve(:,i+1);
                s = s + norm(pj-pi);                
            end
        end
    end
end