classdef triangle < handle
    properties
        n = 0
        pts = zeros(3,2);
        l1
    end
    methods
        function this = triangle()
        end
        function append(this,pt)
            this.n = this.n+1;
            this.pts(this.n,:) = pt;
        end
        function [bool,d] = GetPerpDir2Line(this)
            d = [];
            bool = false;
            A = this.pts(1,:);
            B = this.pts(2,:);
            u = A - B;
            v = -B;
            % are perp?
            v_ortho = [-v(2),v(1)];
            if dot(u,v_ortho) == 0
                bool = true;
            end
            d = this.triple(u,v,u);
            %plot
            this.l1 = this.plotLineAndVector(A,B,d);
        end
        function [bool,d] = ContainOrigin(this,triA,triB)
            A = this.pts(1,:);
            B = this.pts(2,:);
            C = this.pts(3,:);
            CB = B-C;
            CA = A-C;
            OC = -C;
            dCB = this.triple(CA,CB,CB);
            dCA = this.triple(CB,CA,CA);
            %this.plot('green');
            %this.plotLineAndVector(B,C,dCB);
            %this.plotLineAndVector(A,C,dCA);
            if dot(dCB,OC) > 0
                d = dCB;
                this.remove([3,2]);
                triA.remove([3,2]);
                triB.remove([3,2]);
                bool = false;
                delete(this.l1);
                this.l1 = this.plotLineAndVector(B,C,d);
                return
            elseif dot(dCA,OC) > 0
                d = dCA;
                this.remove([3,1]);
                triA.remove([3,1]);
                triB.remove([3,1]);
                bool = false;
                delete(this.l1);
                this.l1 = this.plotLineAndVector(A,C,d);
                return
            end
            
            d = [];
            bool = true;
        end
        function d = triple(this,u,v,z)
            u = [u,0];
            v = [v,0];
            z = [z,0];
            d = cross(cross(u,v),z);
            d = [d(1),d(2)];
        end
        function remove(this,ids)
            this.pts = [this.pts(ids,:);[0,0]];
            this.n = 2;

        end
        function plot(this,color)
            for i = 1:this.n
                line([this.pts(i,1),this.pts(mod(i,this.n)+1,1)],[this.pts(i,2),this.pts(mod(i,this.n)+1,2)],'color',color);
            end
        end
        function lam = getOriginBary(this)
            if this.n == 2
                A = this.pts(1,:);
                B = this.pts(2,:);
                l = norm(B-A);
                lam = [norm(B)/l, norm(A)/l , 0];
            else
                A = [this.pts(1,:),0];
                B = [this.pts(2,:),0];
                C = [this.pts(3,:),0];
                u = B - A;
                v = C - A;
                area = norm(cross(u,v));
                lam1 = norm(cross(-B,-C))/area;
                lam2 = norm(cross(-A,-C))/area;
                lam3 = norm(cross(-A,-B))/area;
                lam = [lam1,lam2,lam3];
                check = lam1 + lam2 + lam3;
            end
        end
        function p = plotBaryPoint(this,lam)
            p = zeros(1,2);
            for i=1:3
                p = p + lam(i)*this.pts(i,:);
            end
            plot(p(1),p(2),'+','color','black','MarkerSize',5);
        end
        function h =plotLineAndVector(this,A,B,d)
            h = line([A(1),B(1)],[A(2),B(2)]);
            c = (A+B)/2;
            n = d/norm(d);
            h(end+1) = quiver(c(1),c(2),n(1),n(2));
        end
    end
end