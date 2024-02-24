classdef Const < FaceGeometry
    properties
        jacobian;
        normal;
    end
    % Constructor
    methods
        function this = Const(fieldNodes,geometryNodes)
            this = this@FaceGeometry(fieldNodes,geometryNodes);
            der = this.getDShapeFunction();
            this.jacobian = norm(cross(der(:,1),der(:,2)));
            this.normal = -(1/(norm(der(:,1))*norm(der(:,1))))* ...
                cross(der(:,1),der(:,2));
        end
    end
    methods
        function [shapeF_Geo,shapeF_For] = getShapeFunction(~,xi,eta)
            shapeF_Geo = zeros(3,1);
            shapeF_Geo(1) = xi;
            shapeF_Geo(2) = eta;
            shapeF_Geo(3) = 1 - xi - eta;
            shapeF_For = [1;1;1];
        end
        function dShapeF = getDShapeFunction(this,~,~)
            dShapeF = zeros(3,2);
            dShapeF(:,1) = this.geometryNodes(1).coord - ...
                this.geometryNodes(3).coord;
            dShapeF(:,2) = this.geometryNodes(2).coord - ...
                this.geometryNodes(3).coord;
        end
                function jacobian = getJacobian(this,~,~)
            jacobian = this.jacobian;
        end
        function normal = getNormal(this,~,~)
            normal = this.normal;
        end
    end
end