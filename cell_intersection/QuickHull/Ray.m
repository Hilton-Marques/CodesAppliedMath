classdef Ray < handle
    properties
        v1
        v2
        dir
    end
    methods
        function this = Ray(v1,v2)
            this.v1 = v1;
            this.v2 = v2;
            u = (v2.coord - v1.coord);
            this.dir = u/norm(u);
        end
        function d = getSquaredDistanceBetweenPoint(this,v)
            c = v - this.v1.coord;
            c_size = dot(c,c);
            a_parrallel = dot(c,this.dir);
            %pythagoras
            d = c_size - a_parrallel*a_parrallel;
        end
        function plot(this)
            pi = this.v1.coord;
            pj = this.v2.coord;
            line([pi(1),pj(1)],[pi(2),pj(2)],[pi(3),pj(3)]);
        end
    end
end