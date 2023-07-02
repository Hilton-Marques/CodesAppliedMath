classdef Cube < handle
    properties
        pts
        n
        centroide
    end
    methods
        function this = Cube(pts)
            this.pts = pts;
            this.n = size(this.pts,1);
            this.centroide = sum(this.pts)/this.n;
        end
        function plot(this,color)
            for i = 1:4
                line([this.pts(i,1),this.pts(mod(i,4)+1,1)],...
                    [this.pts(i,2),this.pts(mod(i,4)+1,2)],[this.pts(i,3),this.pts(mod(i,4)+1,3)],'color',color);
                line([this.pts(i,1),this.pts(i+4,1)],...
                    [this.pts(i,2),this.pts(i+4,2)],[this.pts(i,3),this.pts(i+4,3)],'color',color);
                line([this.pts(i+4,1),this.pts(mod(i,4)+5,1)],...
                    [this.pts(i+4,2),this.pts(mod(i,4)+5,2)],[this.pts(i+4,3),this.pts(mod(i,4)+5,3)],'color',color);
            end
        end
        function transform(this,teta,t)
            rotz = [cos(teta),sin(teta),0;-sin(teta),cos(teta),0;0,0,1];
            this.pts = transpose(rotz*this.pts') + t;
            this.centroide = sum(this.pts)/this.n;
        end
        function v = suportfunction(this,di)
            %pts = this.pts - this.centroide;
            dots = this.pts*di;
            [~, argmax] = max(dots);
            v = this.pts(argmax,:);
        end
        function plot_tri(this)
            conec = [0,1,2;... %up
                1,3,2;...
                
                1,5,7;... %right
                7,3,1;...
                
                4,0,6;... %left
                0,2,6;...
                
                5,1,0;... %front
                0,4,5;...
                
                5,4,6;... %down
                5,6,7;...
                
                7,6,2;...%back
                7,2,3];
            
            conec = conec + 1;
            trisurf(conec,this.pts(:,1), this.pts(:,2), this.pts(:,3),'FaceAlpha',0.2);
        end
    end
end