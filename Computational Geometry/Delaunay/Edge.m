classdef Edge < handle
    properties 
        p1
        p2
        T = Triangle.empty;
    end
    methods
        function this = Edge(p1,p2)
            if nargin > 0
                this.p1 = p1;
                this.p2 = p2;
            end
        end
        function plot(this)
            line([this.p1(1),this.p2(1)],[this.p1(2),this.p2(2)],'Color','y');
        end
    end
end