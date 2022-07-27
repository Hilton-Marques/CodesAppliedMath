classdef Viewer
    %This class is responsible for visualization
    %   Properties:
    %     starNumber, X, Y, Z,markerSize, axisLim, axisLimApply,
    %     minorGridApply, trackApply, highlightApply,
    %   Methods:
    %     Viewer, update
    %   Author: Haihong
    %           Valentine's Day, 2015
    %   inspired by Salman Khan's videos on astronomy at Khan Academy.
    
    properties
        starNumber
        X
        Y
        Z
        markerSize = 20;
        axisLim = 400; % needs a proper value
        axisLimApply = 'off';
        minorGridApply = 'on';
        trackApply = 'on'; % when it's truned of , it is recommended to set axisLimApply to 'on'.
        highlightApply = 'on';
    end
    
    methods
        function obj = Viewer(starNumber,starCluster,t)
            obj.starNumber = starNumber;
            obj.X = 1:obj.starNumber;
            obj.Y = 1:obj.starNumber;
            obj.Z = 1:obj.starNumber;
            obj.update(starCluster,t);
        end
        
        function update(obj,starCluster,t)
            for s = 1:obj.starNumber
                obj.X(s) = starCluster(s).position(1);
                obj.Y(s) = starCluster(s).position(2);
                obj.Z(s) = starCluster(s).position(3);
            end
            axesHandle = gca;
            if strcmp(obj.highlightApply,'on') == true
                scatter3(obj.X(1),obj.Y(1),obj.Z(1),obj.markerSize,'.','m');
                hold(axesHandle,'on');
                scatter3(obj.X(2:obj.starNumber),obj.Y(2:obj.starNumber),obj.Z(2:obj.starNumber),obj.markerSize,'.','b');
            else
                scatter3(obj.X,obj.Y,obj.Z,obj.markerSize,'.','b');
            end
            title(['time = ',num2str(t,'%.4f')]);
            view([40,22]);
            axesHandle.XMinorGrid = obj.minorGridApply;
            axesHandle.YMinorGrid = obj.minorGridApply;
            axesHandle.ZMinorGrid = obj.minorGridApply;
            hold(axesHandle,obj.trackApply);
            xlabel('x');ylabel('y');zlabel('z');
            if strcmp(obj.axisLimApply,'on') == true
                axis(obj.axisLim*[-1,1,-1,1,-1,1]);
            end
        end
    end
end

