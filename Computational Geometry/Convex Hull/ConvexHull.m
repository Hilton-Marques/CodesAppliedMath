classdef ConvexHull < handle
    properties
        points
        N
    end
    methods
        function this = ConvexHull(points)
            this.N = size(points,1);
            this.points = points;
            if this.N == 3
                % reORder points clockwise
                this.reOrder();
            end            
        end
        function plot(this,color)
            points = this.points;
            points(end+1,:) = points(1,:);
            p1 = points(1,:);
            for i = 2:this.N + 1
                plot(p1(1),p1(2),'o','Color','Red');
                p2 = points(i,:);
                line([p1(1),p2(1)],[p1(2),p2(2)],'Color',color);
                p1 = p2; 
            end
        end
        function reOrder(this)
            points = this.points;
            if points(2,2) < points(3,2)
                this.points = this.points([1,3,2],:);
            end
        end
    end
end