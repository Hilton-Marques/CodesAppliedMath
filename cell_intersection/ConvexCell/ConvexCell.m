classdef ConvexCell < handle
    properties
        pts
        n
        centroide
    end
    methods
        function this = ConvexCell(pts)
            this.pts = pts;
            this.n = size(this.pts,1);
            this.centroide = sum(this.pts)/this.n;
        end
        function plot_wireframe(this,color)
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
        function v = supportfcn(this,di)
            dots = (this.pts - this.centroide)*di;
            [~, argmax] = max(dots);
            v = this.pts(argmax,:);
        end
        function show(this,color)
            conec = convhull(this.pts(:,1), this.pts(:,2),this.pts(:,3));
            trisurf(conec,this.pts(:,1), this.pts(:,2), this.pts(:,3),'FaceAlpha',0.6,'FaceColor',color,'EdgeAlpha',1.0);
        end
    end
end