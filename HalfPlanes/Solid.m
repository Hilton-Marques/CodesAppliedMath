classdef Solid < handle
    properties 
        hfs = HalfPlane.empty
        verts;
    end
    methods
        function this = Solid(hfs,verts)
            this.hfs = hfs;
            this.verts = verts;
        end
        function show(this,color)
            if nargin == 1
                color = 'blue';
            end
            for i = 1:length(this.hfs)
                this.hfs(i).show(color);
            end
        end
        function angles = getAngles(this)
            n = length(this.hfs);
            angles = zeros(n,1);
            for i = 1:n
                angles(i) = this.hfs(i).getAngle();
            end            
        end
        function [A,b] = getA(this)
            n = length(this.hfs);
            A = ones(n,4);
            b = zeros(n,1);
            for i = 1:n
                normal = this.hfs(i).getNormal();
                A(i,1:3) = normal;
                b(i) = dot(normal,this.hfs(i).pts(1).coord);
            end
        end
        function pts = getPointsForFecho(this,r)
            pts = zeros(length(this.hfs),3);
            for i = 1:length(this.hfs)
                hfi = this.hfs(i);
                pts(i,:) = hfi.transform();
            end
        end
        function translateHfs(this,center)
            for i = 1:length(this.verts)
                this.verts(i).coord =  this.verts(i).coord - center;
            end
        end
    end
end