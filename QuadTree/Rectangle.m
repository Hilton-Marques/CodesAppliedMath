classdef Rectangle
    properties
        x = [];
        y = [];
        w = [];
        h = [];
    end
    methods
        function this = Rectangle(x,y,w,h) 
            this.x = x;
            this.y = y;
            this.w = w;
            this.h = h;
        end
        function contains = contains(this,point)
            contains = false;
            if( point(1,1) > (this.x - this.w) && ...
                    point(1,1) < (this.x + this.w) && ...
                    point(1,2) > (this.y - this.h) && ...
                    point(1,2) < (this.y + this.h))
                contains = true;
            end
            
        end
    end
end