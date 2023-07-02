classdef Quadrado < handle
    properties
        points
        heds = zeros(4,4);
        id;
    end
    methods
        function this = Quadrado(points,id)
            if nargin > 0
                this.points = points;
                this.createHeds();
                this.id = id;
            end
        end
        function createHeds(this)
            for i = 1:4
                p1 = this.points(i,:);
                p2 = this.points(mod(i,4)+1,:);
                hed = [p1,p2];
                this.heds(i,:) = hed;
            end
        end
        function bool = isInside(this,pt)
            bool = false;
            if (pt(1) >= this.points(1,1) && pt(1) <= this.points(2,1)) && ...
                    (pt(2) >= this.points(1,2) && pt(2) <= this.points(3,2))
                bool = true;
            end
        end
        function plot(this)
            for i = 1:4
                hed = this.heds(i,:);
                line([hed(1),hed(3)],[hed(2),hed(4)],'Color','blue');
            end
        end
    end
end