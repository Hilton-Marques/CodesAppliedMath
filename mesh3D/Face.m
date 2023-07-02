classdef Face < handle
    properties
        pts;
        nChild = 4;
        child = cell(4,1);
        level;
        isDivided = false;
        key = '0';
    end
    methods
        function this = Face(points,level,key)
            if nargin > 0
                this.pts = points;
                this.level = level;
                this.key = key;
            end
        end
        function subdivide(this)
            this.isDivided = true;
            p1 = this.pts(1,:);
            p2 = this.pts(2,:);
            p3 = this.pts(3,:);
            p4 = this.pts(4,:);
            pmed1 = this.mean(p1,p2);
            pmed2 = this.mean(p2,p3);
            pmed3 = this.mean(p3,p4);
            pmed4 = this.mean(p4,p1);
            pcg = this.mean(p1,p3);
            nextLevel = this.level + 1;
            this.child{1} = Face([p1;pmed1;pcg;pmed4],nextLevel,strcat(this.key,'0'));
            this.child{2} = Face([pmed1;p2;pmed2;pcg],nextLevel,strcat(this.key,'1'));
            this.child{3} = Face([pcg;pmed2;p3;pmed3],nextLevel,strcat(this.key,'2'));
            this.child{4} = Face([pmed4;pcg;pmed3;p4],nextLevel,strcat(this.key,'3'));
        end
        function out = mean(~,p1,p2)
            out = (p1+p2)/2;
        end
        function show(this,off)
            p1 = this.pts(1,:);
            p2 = this.pts(2,:);
            p3 = this.pts(3,:);
            p4 = this.pts(4,:);
            n = cross((p2 - p1),(p4-p1));
            n = n/norm(n);
            p1 = p1 + off*n;
            p2 = p2 + off*n;
            p3 = p3 + off*n;
            p4 = p4 + off*n;            
            line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)]);
            line([p2(1),p3(1)],[p2(2),p3(2)],[p2(3),p3(3)]);
            line([p3(1),p4(1)],[p3(2),p4(2)],[p3(3),p4(3)]);
            line([p4(1),p1(1)],[p4(2),p1(2)],[p4(3),p1(3)]);
        end
    end
end