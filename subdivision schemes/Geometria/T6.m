classdef T6 < FaceGeometry
    % Constructor
    methods
        function this = T6(fieldNodes,geometryNodes)
            this = this@FaceGeometry(fieldNodes,geometryNodes);
        end
    end
    methods
        function [shapeF_Geo,shapeF_For] = getShapeFunction(~,xi,eta)
            shapeF_Geo = zeros(6,1);
            zeta = 1 - xi - eta;
            shapeF_Geo(1) = (2*xi-1)*xi;
            shapeF_Geo(2) = (2*eta - 1)*eta;
            shapeF_Geo(3) = (2*zeta - 1)*zeta;
            shapeF_Geo(4) = 4*xi*eta;
            shapeF_Geo(5) = 4*eta*zeta;
            shapeF_Geo(6) = 4*zeta*xi;
            shapeF_For = shapeF_Geo;
        end
        function dShapeF = getDShapeFunction(this,xi,eta)
            dShapeF = zeros(3,2);
            zeta = 1 - xi - eta;
            dShapeF(:,1) = 4*eta*this.geometryNodes(4).coord + ...
                2*this.geometryNodes(1).coord*xi + ...
                this.geometryNodes(1).coord*((-1)+2*xi)+...
                +4*this.geometryNodes(6).coord*zeta;
            dShapeF(:,2) = 2*eta*this.geometryNodes(2).coord + ...
                ((-1)+2*eta)*this.geometryNodes(2).coord + ...
                4*this.geometryNodes(4).coord*xi + ...
                4*this.geometryNodes(5).coord*zeta;
        end
        function jacobian = getJacobian(this,xi,eta)
            der = this.getDShapeFunction(xi,eta);
            jacobian = norm(cross(der(:,1),der(:,2)));
        end
        function normal = getNormal(this,xi,eta)
            der = this.getDShapeFunction(xi,eta);
            normal = -(1/(norm(der(:,1))*norm(der(:,1))))* ...
                cross(der(:,1),der(:,2));
        end
    end
end