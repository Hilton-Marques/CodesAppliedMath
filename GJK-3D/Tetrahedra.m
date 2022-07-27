classdef Tetrahedra < handle
    properties
        n = 0
        pts = zeros(4,3);
        tri_image = []
    end
    methods
        function this = Tetrahedra()
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
            u_ortho = this.getPtOnPerpPlane(u);
            if dot(v,u_ortho) == 0
                bool = true;
            end
            d = this.triple(u,v,u);
            %plot
            this.plotLineAndVector(A,B,d);
        end
        function [bool,d] = GetPerpDir2Tri(this)
            d = [];
            bool = false;
            A = this.pts(1,:);
            B = this.pts(2,:);
            C = this.pts(3,:);
            u = C - A;
            v = C - B;
            d = cross(u,v);
            dotnc = dot(d,-C);
            if dotnc == 0
                bool = true;
            elseif dotnc < 0
                d = -d;
            end
            this.plotTriAndVector(A,B,C,d);
        end
        function x = getPtOnPerpPlane(this,line)
            z = line/norm(line);
            uTemp = [z(3);z(1);-z(2)];
            uTemp = uTemp/norm(uTemp);
            if (dot(z,uTemp) == 1)
                uTemp = uTemp([2,1,3]);
            end
            x = cross(uTemp,z);
        end
        function d = triple(this,u,v,z)
            d = cross(cross(u,v),z);
        end
        function plotLineAndVector(this,A,B,d)
            len = norm(B-A);
            line([A(1),B(1)],[A(2),B(2)],[A(3),B(3)]);
            c = (A+B)/2;
            n = (len/3)*d/norm(d);
            quiver3(c(1),c(2),c(3),n(1),n(2),n(3));
        end
        function plotTriAndVector(this,A,B,C,d)
            pts = [A;B;C];
            area = norm(cross(B-A,C-A));
            d = (sqrt(area)/3)* d/norm(d);
            this.tri_image = trisurf([1,2,3],pts(:,1),pts(:,2),pts(:,3),'FaceAlpha',0.5);
            centroide = sum(pts,1)/3;
            this.tri_image(end+1) = quiver3(centroide(1),centroide(2),centroide(3),d(1),d(2),d(3),'color','black');
        end
        function lam = getOriginBary(this)
            if this.n == 2
                A = this.pts(1,:);
                B = this.pts(2,:);
                l = norm(B-A);
                lam = [norm(B)/l, norm(A)/l , 0,0];
            elseif this.n == 3
                A = this.pts(1,:);
                B = this.pts(2,:);
                C = this.pts(3,:);
                u = B - A;
                v = C - A;
                area = norm(cross(u,v));
                lam1 = norm(cross(-B,-C))/area;
                lam2 = norm(cross(-A,-C))/area;
                lam3 = norm(cross(-A,-B))/area;
                lam = [lam1,lam2,lam3,0];
            else
                A = this.pts(1,:);
                B = this.pts(2,:);
                C = this.pts(3,:);
                D = this.pts(4,:);
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
            A = this.pts(1,:);
            B = this.pts(2,:);
            C = this.pts(3,:);
            D = this.pts(4,:);
            OD = -D;
            dDAB = this.tripleTri([4,1,2],3);
            dDBC = this.tripleTri([4,2,3],1);
            dDCA = this.tripleTri([4,3,1],2);
            figure 
            hold on
            view(30,30);
            lighting gouraud;
            h = this.plotTetra();
            quiver3(D(1),D(2),D(3),OD(1),OD(2),OD(3));
            plot3(0,0,0,'o','MarkerFaceColor','cyan','MarkerSize',5);
            p = OD / norm(OD);
            dDAB = dDAB/norm(dDAB);
            dDBC = dDBC/norm(dDBC);
            dDCA = dDCA/norm(dDCA);
            values = [dDAB;dDBC;dDCA]*p'
            [max_value,argmax] = max([dDAB;dDBC;dDCA]*p');
            
            if (max_value > 0)
            %if dot(dDAB,OD) > 0
            if argmax == 1
                d = dDAB;
                %dot_prod = (dot(OD,d)/(norm(OD)*norm(d)))
                this.remove([4,1,2]);
                triA.remove([4,1,2]);
                triB.remove([4,1,2]);
                delete(this.tri_image);
                bool = false;
                this.plotTriAndVector(D,A,B,d);
                delete(h);
                hold off
                return
            %elseif dot(dDBC,OD) > 0
            elseif argmax == 2
                d = dDBC;
                dot_prod = (dot(OD,d)/(norm(OD)*norm(d)))
                this.remove([4,2,3]);
                triA.remove([4,2,3]);
                triB.remove([4,2,3]);
                delete(this.tri_image);
                bool = false;
                this.plotTriAndVector(D,B,C,d);
                delete(h);
                hold off
                return
            %elseif dot(dDCA,OD) > 0
            elseif argmax == 3
                d = dDCA;
                dot_prod = (dot(OD,d)/(norm(OD)*norm(d)))
                this.remove([4,3,1]);
                triA.remove([4,3,1]);
                triB.remove([4,3,1]);
                delete(this.tri_image);
                bool = false;
                this.plotTriAndVector(D,A,C,d);
                delete(h);
                hold off
                return
            end
            end
            d = [];
            bool = true;
        end
        function p = plotBaryPoint(this,lam)
            p = zeros(1,3);
            for i=1:4
                p = p + lam(i)*this.pts(i,:);
            end
            %plot(p(1),p(2),'+');
        end
        function remove(this,ids)
            pts = this.pts(ids,:);
            this.pts = zeros(4,3);
            this.n = 3;
            this.pts(1:3,:) = pts;
        end
        function d = tripleTri(this,ids,i)
            A = this.pts(ids(1),:);
            B = this.pts(ids(2),:);
            C = this.pts(ids(3),:);
            O = this.pts(i,:);
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
            A = this.pts(1,:);
            B = this.pts(2,:);
            C = this.pts(3,:);
            D = this.pts(4,:);
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
        function pt = pertubate(this,pt,coords,normal,lam)
            for i = 1:4
                A = this.pts(i,:);
                for j = (i+1):4
                    B = this.pts(j,:);
                    if (norm(A - B) == 0)
                        lam(i) = lam(i) + lam(j);
                        lam(j) = 0;
                        break;
                    end
                end
            end
            value = 0;
            id = -1;
            for i = 1:4
                if (lam(i) > value)
                    value = lam(i);
                    id = i;
                end
            end
            A = this.pts(id,:);
            plot3(A(1),A(2),A(3),'o','MarkerSize',10);
            %find A on coords
            for i = 8:-1:1
                if norm((coords(i,:) - A)) == 0
                    id = i;
                    break;
                end
            end
            opp = coords(8 - id + 1,:);
            if id < 5
                normal = -normal;
            end
            dir = opp + 0.1*normal - pt;
            dir = dir/norm(dir);
            pt = pt + 0.1*dir;
        end
        function h = plotTetra(this)
            k = [1,2,3;1,2,4;2,3,4;3,1,4];
            h = trisurf(k,this.pts(:,1),this.pts(:,2),this.pts(:,3),'FaceAlpha',0.5);
        end
        
    end
end