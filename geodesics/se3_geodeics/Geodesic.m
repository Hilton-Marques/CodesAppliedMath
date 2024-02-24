classdef Geodesic < Scene
    properties

    end
    methods 
        function this = Geodesic(A,B)
            this@Scene('geodesic_se3.gif','SE3');
            v = A.logMap(A.m_data, B.m_data);
            n = 10;
            t = linspace(0,1,n);
            pts = A.geodesic(v,A.m_data,t);
            coords = pts(:,1:3,4)';
            this.setMargin(2);
            this.setBB(coords);
            view(-45,25);
            this.drawRobot(squeeze(pts(1,:,:)), color=this.m_blue);
            this.drawBoard();
            %this.plotText(pts(1,1:3,4) - 2* [1,0,0], '$\mathcal{X}$');
            %this.get();
            %this.exportFrame('A');
            %this.drawRobot(squeeze(pts(end,:,:)), this.m_red, 0.8);
            %this.get();
            %this.plotText(pts(end,1:3,4) - 2* [1,0, 0], '$\mathcal{Y}$');
            %this.exportFrame('B');
            for i = 2:n-1
                color = this.m_blue + (i/10) * (this.m_red - this.m_blue);
                h = this.drawRobot(squeeze(pts(i,:,:)),color=color);
                %pause(0.1);
                this.get();
                %delete(h);
            end
            this.exportFrame();
            this.save(true)
        end
    end
end