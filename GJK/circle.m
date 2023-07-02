classdef circle < handle
    properties
        pts
        n
        centroide
    end
    methods
        function this = circle(r)
            t = linspace(0,2*pi,30);
            this.pts = r*[cos(t)',sin(t)'];
            this.n = size(this.pts,1);
            this.centroide = sum(this.pts)/this.n;
        end
        function plot(this,color)
            for i = 1:this.n
                line([this.pts(i,1),this.pts(mod(i,this.n)+1,1)],[this.pts(i,2),this.pts(mod(i,this.n)+1,2)],'color',color);
            end
        end
        function v = suportfunction(this,di)
            pts = this.pts - this.centroide;
            dots = pts*di;
            [~, argmax] = max(dots);
            v = this.pts(argmax,:);
        end
        function transform(this,teta,t)
            R = [cos(teta),-sin(teta);sin(teta),cos(teta)];
            this.pts = transpose(R*this.pts') + t;
            this.centroide = sum(this.pts)/this.n;
        end
    end
end