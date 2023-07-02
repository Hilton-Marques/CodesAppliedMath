classdef Element < handle
    properties
        index = [];
        points = [];
        face = [];
        q = [];
        id = 0;
        area = [];
        n = [];
    end
    methods
        function this = Element(index,points,face,i)
            if nargin > 0
                this.index = index;
                this.points = points(index,:); 
                this.face = face;
                this.q = face.q;
                this.id = i;    
                this.findArea();   
            end
        end
        function findArea(this)
            u = this.points(2,:) - this.points(1,:);
            v = this.points(3,:) - this.points(1,:);
            biV = cross(u,v);
            this.n = biV/norm(biV);
            this.area = norm(biV)*0.5;
        end
    end
end