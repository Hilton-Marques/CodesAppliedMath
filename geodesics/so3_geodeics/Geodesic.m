classdef Geodesic < Scene
    properties

    end
    methods 
        function this = Geodesic(A,B)
            this@Scene('geodesic_so3.gif','SO3');
            v = A.logMap(A.m_data, B.m_data);
            n = 10;
            t = linspace(0,1,n);
            pts = A.geodesic(v,A.m_data,t);
            p0 = zeros(3,1);
            p1 = n*[1; 0; 0];
            coords = p0 + linspace(0,1,n).*(p1 - p0);
            %coords = pts(:,1:3,4)';
            %this.setMargin(2);
            %this.setBB(coords);
            view(-45,25);
            T = squeeze(pts(1,:,:));
            T(1:3,4) = p0;
            T(4, 4) = 1.0;
            this.drawAxis(T, this.m_blue);
            %this.drawRobot(squeeze(pts(1,:,:)), this.m_blue, 0.8);
            %this.plotText(p0, '$\mathcal{X}$');
            %this.get();
            %this.exportFrame('A');
            %this.drawRobot(squeeze(pts(end,:,:)), this.m_red, 0.8);
            %this.get();
            %this.plotText(p1, '$\mathcal{Y}$');
            %this.exportFrame('B');
            for i = 2:n
                color = this.m_blue + (i/10) * (this.m_red - this.m_blue);
                T = squeeze(pts(i,:,:));
                T(1:3, 4) = coords(1:3,i);
                T(4,4) = 1.0;
                %h = this.drawRobot(squeeze(pts(i,:,:)),color,0.8);
                h = this.drawAxis(T, color);
                %pause(0.1);
                this.get();
                %delete(h);
            end
            this.exportFrame();
            this.save(true)
        end
    end
end