classdef FaceGeometry < handle
    properties
        fieldNodes;
        geometryNodes;
    end
    methods
        % Method used in BEM
        function this = FaceGeometry(fieldNodes,geometryNodes)
            this.fieldNodes = fieldNodes;
            this.geometryNodes = geometryNodes;
        end
        %These functions could be obtained here, but as these quatities
        %are constant for Const and T3, I put on her classes.
        
        %         function jacobian = getJacobian(this,xi,eta)
        %             der = this.getDShapeFunction(xi,eta);
        %             jacobian = norm(cross(der(:,1),der(:,2)));
        %         end
        %         function normal = getNormal(this,xi,eta)
        %             der = this.getDShapeFunction(xi,eta);
        %             normal = -(1/(norm(der(:,1))*norm(der(:,1))))* ...
        %                 cross(der(:,1),der(:,2));
        %         end
        
        % method used in FMM
        function yc = getYc(this)
            yc = zeros(3,1);
            n = 3;
            for i = 1:n
                yc = yc + this.geometryNodes(i).coord;
            end
            yc = yc/n;
        end
        function drawCircumSphere(this,fac)
            u = this.geometryNodes(1).coord - this.geometryNodes(2).coord ;
            v = this.geometryNodes(1).coord - this.geometryNodes(3).coord ;
            n = cross(u,v);
            n = n/norm(n);
            ortho_u = cross(n,u);
            ortho_v = cross(n,v);
            a = u;
            b = v;
            ac = ortho_u;
            bc = ortho_v;
            rc = (dot(b,b)*ac - dot(a,a)*bc)/(2*dot(ac,b));
            radius = fac*norm(rc);
            center = this.geometryNodes(1).coord - rc;
            [X,Y,Z]  = sphere();
            X2 = radius * X + center(1);
            Y2 = radius * Y + center(2);
            Z2 = radius * Z + center(3);
            surf(X2,Y2,Z2,'FaceAlpha',0.3,'FaceColor','red','EdgeAlpha',0.1);
        end
        
        function radius = getRadius(this,fac)
            if nargin == 1
                fac = 1;
            end
            u = this.geometryNodes(1).coord - this.geometryNodes(2).coord ;
            v = this.geometryNodes(1).coord - this.geometryNodes(3).coord ;
            z = this.geometryNodes(2).coord - this.geometryNodes(3).coord ;
            L1 = dot(u,u);
            L2 = dot(v,v);
            L3 = dot(z,z);
            if (L1 > L2)
                if L1 > L3
                    radius = L1;
                else
                    radius = L3;
                end
            else
                if L2 > L3
                    radius = L2;
                else
                    radius = L3;
                end
            end
            radius = fac^2  * radius;
        end
        function plotSpheres(this,fac)
            radius = sqrt(this.getRadius(fac));
            [X,Y,Z]  = sphere();
            for i = 1:3
                trans = this.geometryNodes(i).coord;
                X2 = radius * X + trans(1);
                Y2 = radius * Y + trans(2);
                Z2 = radius * Z + trans(3);
                surf(X2,Y2,Z2,'FaceAlpha',0.3,'FaceColor','cyan','EdgeAlpha',0.1);
            end
            
        end
    end
    methods (Abstract)
        [shapeF_Geo,shapeF_For] = getShapeFunction(~,xi,eta);
        dShapeF = getDShapeFunction(this,xi,eta);
    end
end