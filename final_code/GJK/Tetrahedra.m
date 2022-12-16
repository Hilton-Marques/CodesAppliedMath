classdef Tetrahedra < handle
    properties
        m_n = 0
        m_pts = zeros(4,3)
        m_tol = 1.9e-7
    end
    methods
        function this = Tetrahedra()
        end
        function append(this,pt)
            this.m_n = this.m_n+1;
            this.m_pts(this.m_n,:) = pt;
        end
        function [bool,d] = GetPerpDir2Line(this)
            d = [];
            bool = false;

            A = this.m_pts(1,:);
            B = this.m_pts(2,:);
            u = A - B;
            v = -B;
            % are perp?
            u_ortho = this.getPtOnPerpPlane(u);
            if abs(dot(v,u_ortho)) < this.m_tol
                bool = true;
            end
            d = this.triple(u,v,u);
        end
        function [bool,d] = GetPerpDir2Tri(this)
            d = [];
            bool = false;
            A = this.m_pts(1,:);
            B = this.m_pts(2,:);
            C = this.m_pts(3,:);
            u = C - A;
            v = C - B;
            d = cross(u,v);
            dotnc = dot(d,-C);
            if dotnc == 0
                bool = true;
            elseif dotnc < 0
                d = -d;
            end
        end
        function p = getBaryPoint(this,lam)
            p = zeros(1,3);
            for i=1:4
                p = p + lam(i)*this.m_pts(i,:);
            end
        end
        function lam = getOriginBary(this)
            if this.m_n == 2
                A = this.m_pts(1,:);
                B = this.m_pts(2,:);
                l = norm(B-A);
                lam = [norm(B)/l, norm(A)/l , 0,0];
            elseif this.m_n == 3
                A = this.m_pts(1,:);
                B = this.m_pts(2,:);
                C = this.m_pts(3,:);
                u = B - A;
                v = C - A;
                area = norm(cross(u,v));
                lam1 = norm(cross(-B,-C))/area;
                lam2 = norm(cross(-A,-C))/area;
                lam3 = norm(cross(-A,-B))/area;
                lam = [lam1,lam2,lam3,0];
            else
                A = this.m_pts(1,:);
                B = this.m_pts(2,:);
                C = this.m_pts(3,:);
                D = this.m_pts(4,:);
                n1 = cross(B-A,C-A);
                n2 = cross(D-A,B-A);
                n3 = cross(D-A,D-C);
                n4 = cross(D-B,D-C);
                volume = norm(dot(n1,D-A));
                lam1 = norm(dot(n4,-D))/volume;
                lam2 = norm(dot(n3,-D))/volume;
                lam3 = norm(dot(n2,-D))/volume;
                lam4 = norm(dot(n1,-A))/volume;
                lam = [lam1,lam2,lam3,lam4];
            end
        end
        function [bool,d] = ContainOrigin(this,triA,triB)
            A = this.m_pts(1,:);
            B = this.m_pts(2,:);
            C = this.m_pts(3,:);
            D = this.m_pts(4,:);
            OD = -D;
            dDAB = this.tripleTri([4,1,2],3);
            dDBC = this.tripleTri([4,2,3],1);
            dDCA = this.tripleTri([4,3,1],2);
            p = OD / norm(OD);
            dDAB = dDAB/norm(dDAB);
            dDBC = dDBC/norm(dDBC);
            dDCA = dDCA/norm(dDCA);
            values = [dDAB;dDBC;dDCA]*p';
            [max_value,argmax] = max([dDAB;dDBC;dDCA]*p');
            if (max_value >= -this.m_tol)
                %if dot(dDAB,OD) > 0
                if argmax == 1
                    d = dDAB;
                    %dot_prod = (dot(OD,d)/(norm(OD)*norm(d)))
                    this.remove([4,1,2]);
                    triA.remove([4,1,2]);
                    triB.remove([4,1,2]);
                    bool = false;
                    return
                elseif argmax == 2
                    d = dDBC;
                    bool = false;

                    dot_prod = (dot(OD,d)/(norm(OD)*norm(d)));
                    this.remove([4,2,3]);
                    triA.remove([4,2,3]);
                    triB.remove([4,2,3]);
                    return
                    %elseif dot(dDCA,OD) > 0
                elseif argmax == 3
                    d = dDCA;
                    bool = false;

                    dot_prod = (dot(OD,d)/(norm(OD)*norm(d)));
                    this.remove([4,3,1]);
                    triA.remove([4,3,1]);
                    triB.remove([4,3,1]);                   
                    return
                end
            end
            d = [];
            bool = true;
        end
        function remove(this,ids)
            pts = this.m_pts(ids,:);
            this.m_pts = zeros(4,3);
            this.m_n = 3;
            this.m_pts(1:3,:) = pts;
        end
        function d = tripleTri(this,ids,i)
            A = this.m_pts(ids(1),:);
            B = this.m_pts(ids(2),:);
            C = this.m_pts(ids(3),:);
            O = this.m_pts(i,:);
            u = C - A;
            v = C - B;
            d = cross(u,v);
            AO = O - A;
            if dot(d,AO) > 0
                d = -d;
            end
        end
        function out = checkDegenerate(this)
            out = [0,0,0];
            count = 1;
            A = this.m_pts(1,:);
            B = this.m_pts(2,:);
            C = this.m_pts(3,:);
            D = this.m_pts(4,:);
            c_pts = A;
            if ~(norm(B - A) == 0)
                count = count + 1;
                c_pts(end+1,:) = B;
            end
            if ((norm(C - A) ~= 0) && norm(C-B) ~= 0)
                count = count + 1;
                c_pts(end+1,:) = C;
            end
            if ((norm(D - A) ~= 0) && (norm(D - B) ~= 0) && (norm(D-C) ~= 0))
                count = count + 1;
                c_pts(end+1,:) = D;
            end
            if (count == 3)
                u = c_pts(2,:) - c_pts(1,:);
                v = c_pts(3,:) - c_pts(1,:);
                n = cross(u,v);
                n = n/norm(n);
                out = n;
            end
        end
        function h = plotTetra(this)
            k = [1,2,3;1,2,4;2,3,4;3,1,4];
            h = trisurf(k,this.m_pts(:,1),this.m_pts(:,2),this.m_pts(:,3),'FaceAlpha',0.5);
        end
        function bool = isInside(this)
            pt = [0,0,0];
            bool = true;
            k = convhull(this.m_pts(:,1),this.m_pts(:,2),this.m_pts(:,3));
            n = size(k,1);
            for i = 1:n
                p0 = this.m_pts(k(i,1),:);
                p1 = this.m_pts(k(i,2),:);
                p2 = this.m_pts(k(i,3),:);
                n = -cross(p2 - p0, p1 - p0);
                u = pt - p0;
                if dot(n,u) > 0
                    bool = false;
                    return;
                end
            end
        end
    end
    methods (Static)
        function x = getPtOnPerpPlane(line)
            z = line/norm(line);
            uTemp = [z(3);z(1);-z(2)];
            uTemp = uTemp/norm(uTemp);
            if (dot(z,uTemp) == 1)
                uTemp = uTemp([2,1,3]);
            end
            x = cross(uTemp,z);
        end
        function d = triple(u,v,z)
            d = cross(cross(u,v),z);
        end
    end
end