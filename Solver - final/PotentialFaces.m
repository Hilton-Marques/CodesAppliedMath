
%
%% Class definition
classdef PotentialFaces < handle
    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        k = 0;        % Conductivty
        len = 0;           % member length
        q = [];             % Boundary Condition
        cord = [];      % Coordinate of each element
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function element = PotentialFaces(k,len,q,cord)
            if (nargin > 0)
                element.k = k;
                element.len = len;
                element.q = q;
                element.cord = cord; 
            
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of a CrossMember object.
        function memb = clean(memb)
            memb.len = 0;
            memb.q = 0;
            memb.k = 0;
            memb.cord = 0; 
        end
    end
end