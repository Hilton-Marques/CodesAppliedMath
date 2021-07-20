classdef TriPolygon < handle
    properties
        polygon = [];
        N = 0;
        index = [];
        L = [];
        R = [];
        triIndex = [];
    end
    methods
        function this = TriPolygon(points,N)
            this.polygon = points;
            this.N = N;
            this.index = 1:N;
            this.L = [N,1:N-1];
            this.R = [2:N,1];
        end
        function startTri(this)
            for i = 1:this.N
                if (this.isAEar(i))
                    this.triIndex(end+1,:) = [this.L(i),i,this.R(i)];
                    l = this.L(i);
                    r = this.R(i);
                    this.L(r) = l;
                    this.R(l) = r;
                end
            end
        end
        function out = isAEar(this,i)
            j = this.L(i);
            k = this.R(i);
            p1 = this.polygon(j,:);
            p2 = this.polygon(i,:);
            p3 = this.polygon(k,:);
            triangle = [p1;p2;p3];
            if (this.signedArea(p1,p2,p3) < 0)
                out = false;
                return
            end
            for l = 1:this.N
                if ((l ~= i) && (l~=j) && (l~=k))
                    p = this.polygon(l,:);
                    if (this.isPointInTriangle(p,triangle))
                        out = false;
                        return;
                    end
                end
            end
            out = true;
        end
        function out = isPointInTriangle(this,c,triangle)
            index = [1,2,3,1];
            for i = 1:3
                a = triangle(index(i),:);
                b = triangle(index(i+1),:);
                if (this.signedArea(a,b,c) < 0)
                    out = false;
                    return
                end
            end
            out = true;
        end
        function out = signedArea(~,a,b,c)
            u = (b - a)/vecnorm(b-a);
            v = (c - a)/vecnorm(c-a);
            out = det([u;v]);
            area = abs(out);
            if area < 0.02
                out = 0;
            end
        end
    end
end