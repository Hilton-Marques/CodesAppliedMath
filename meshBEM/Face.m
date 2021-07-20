classdef Face < handle
    properties
        hedInc;
        pointsInc;
    end
    methods
        function this = Face(hedInc,heds)
            if nargin > 0
                this.hedInc = hedInc;
                this.getIndcPoints(heds)
            end
        end
        function getIndcPoints(this,heds)
            hed1 = heds(this.hedInc(1));
            p1 = hed1.inc(1);
            p2 = hed1.inc(2);
            hed2 = heds(this.hedInc(2));
            ptemp = hed2.inc(1);
            if ptemp == p2
                p3 = hed2.inc(2);
            else
                p3 = hed2.inc(1);
            end
            this.pointsInc = [p1,p2,p3];
        end
        function plot(this,points)
            p = points(this.pointsInc);
            line([p(1).coord(1),p(2).coord(1)],[p(1).coord(2),p(2).coord(2)],...
                [p(1).coord(3),p(2).coord(3)]);
            line([p(2).coord(1),p(3).coord(1)],[p(2).coord(2),p(3).coord(2)],...
                [p(2).coord(3),p(3).coord(3)]);
            line([p(3).coord(1),p(1).coord(1)],[p(3).coord(2),p(1).coord(2)],...
                [p(3).coord(3),p(1).coord(3)]);
            p = [p(1).coord';p(2).coord';p(3).coord'];
            plot3(p(:,1), p(:,2), p(:,3),'o');
        end
    end
end