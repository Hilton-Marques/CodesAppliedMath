classdef Geodesic < Scene
    properties

    end
    methods 
        function this = Geodesic(A,B)
            this@Scene('geodesic_se2.gif','SE2');
            v = A.logMap(A.m_data, B.m_data);
            n = 10;
            t = linspace(0,1,n);
            pts = A.geodesic(v,A.m_data,t);
            coords = pts(:,1:2,3)';
            this.setMargin(3);
            coords(3,:) = 0;
            this.setBB(coords);
            this.drawRobot(squeeze(pts(1,:,:)), this.m_blue);
            %this.plotText2D(pts(1,1:2,3) - 0.5* [1,0], '$\mathcal{X}$');
            %this.get();
            %this.exportFrame('A');
            this.drawRobot(squeeze(pts(end,:,:)), this.m_red);
            %this.get();
            tf = B.m_data(1:2, end);
            %this.plotText2D(tf + [0;0.25], '$\mathcal{Y}$');
            %this.exportFrame('B');

            for i = 2:n-1
                color = this.m_blue + (i/10) * (this.m_red - this.m_blue);
                h = this.drawRobot(squeeze(pts(i,:,:)),color);
                %pause(0.1);
                this.get();
                %delete(h);
            end
            this.exportFrame();
            this.save(true)
        end
    end
end