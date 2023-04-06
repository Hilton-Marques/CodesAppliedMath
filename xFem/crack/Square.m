classdef Square < handle
    properties
        pts
        heatmap
        inter = 0
    end
    methods
        function this = Square(pts,heatmap)
            if nargin > 0
                this.pts = pts;
                this.heatmap = heatmap;
            end
        end
        function bool = crackIntersection(this,crack)
            bool = false;
            for i = 1:4
                xi = this.pts(i,:);
                xj = this.pts(mod(i,4)+1,:);
                if this.lineline([xi,xj],crack)
                    bool = true;
                    this.inter = this.inter + 1;
                end
            end
        end
        function show(this)
            pts = this.pts;
            pts(:,3) = 0;
%             X = zeros(2,2);
%             X(1,1) = this.pts(1,1);
%             X(1,2) = this.pts(2,1);
%             X(2,1) = this.pts(4,1);
%             X(2,2) = this.pts(3,1);
%             Y = zeros(2,2);
%             Y(1,1) = this.pts(1,2);
%             Y(1,2) = this.pts(2,2);
%             Y(2,1) = this.pts(4,2);
%             Y(2,2) = this.pts(3,2);
%             Z = reshape(this.heatmap,2,2);
            %pcolor(X,Y,Z);
            trisurf([1,2,3;3,4,1],pts(:,1),pts(:,2),pts(:,3),'FaceAlpha',1.0,'FaceColor','white','EdgeAlpha',0.0);
        end
        function bool = lineline(this,hed,line)
            bool = false;
            C = line(1:2);
            D = line(3:4);
            A = hed(1:2);
            B = hed(3:4);
            if this.orient(C,D,A) > 0 && this.orient(C,D,B) > 0
                return
            elseif this.orient(C,D,A) < 0 && this.orient(C,D,B) < 0
                return
            elseif this.orient(A,B,C) > 0 && this.orient(A,B,D) > 0
                return
            elseif this.orient(A,B,C) < 0 && this.orient(A,B,D) < 0
                return
            end
            bool = true;
            return
        end
    end
    methods (Static)
        
        function out = orient(A,B,C)
            out = det([B-A;C-A]);
        end
    end
end