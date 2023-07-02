classdef Triangle < handle
    properties
        p1;
        p2;
        p3;
    end
    methods
        function this = Triangle(p1,p2,p3)
            this.p1 = p1;
            this.p2 = p2;
            this.p3 = p3;
        end
        function show(this)
            fill3([this.p1(1),this.p2(1), this.p3(1)], ...
                [this.p1(2),this.p2(2), this.p3(2)], ...
                [this.p1(3),this.p2(3), this.p3(3)], [0,0,1], 'faceAlpha',0.5); 
        end
        function centroide = getCentroide(this)
            centroide = (this.p1 + this.p2 + this.p3)/3 ;
        end
    end
end