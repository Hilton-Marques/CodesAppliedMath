classdef Element < handle
    properties
        k = 1; %% conductivity
        area;  %% area
        q = []; %% boundary
        nNos;  %% nodes per element
        face = []; %% face where the element belong
        pos = []; %% coordinates
        id = 1;   %% identification
        n;  %% normal vector
    end
    methods
        function this = Element(face,nNos,index,points,id)
            if nargin > 0
                this.face = face;
                this.nNos = nNos;
                this.pos = points(index,:);
                this.q = face.q;
                this.id = id;
                this.findArea;
            end
            
        end
        function findArea(this)
            u = this.pos(2,:) - this.pos(1,:);
            v = this.pos(3,:) - this.pos(1,:);
            biV = cross(u,v);
            this.n = biV'/norm(biV);
            this.area = norm(biV)*0.5;
        end
        function plot(this)
            p1 = this.pos(1,:);
            p2 = this.pos(2,:);
            p3 = this.pos(3,:);
            line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'Color','black');
            line([p2(1),p3(1)],[p2(2),p3(2)],[p2(3),p3(3)],'Color','black');
            line([p3(1),p1(1)],[p3(2),p1(2)],[p3(3),p1(3)],'Color','black');
        end
    end
end