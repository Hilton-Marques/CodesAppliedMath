
%
%% Class definition
classdef PotentialElement < handle
    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        k = 0;             % Conductivty
        len = 0;           % member length
        q = [];             % Boundary Condition
        nnos = 0;          % Nodes per element
        face = [];
        cord = [];      % Coordinate of each element
        indnodes = [];
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------

        function init(element)
            element.k = element.face.k;
            element.nnos = numel(element.indnodes);
            element.q = element.face.q;
            element.len = norm([element.cord(3)-element.cord(1); ...
                         element.cord(4) - element.cord(2)]);
        end
    end
    
    %% Public methods

end