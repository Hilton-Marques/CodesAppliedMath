
%
%% Class definition
classdef PotentialElement < handle
    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        k = 0;             % Conductivty
        len = 0;           % member length
        q = [];             % Boundary Condition
        nnos = 0;          % Nodes per element
        face = 0;
        cord = [];      % Coordinate of each element
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function element = PotentialElement(k,len,q,nnos,face,cord)
            if (nargin > 0)
                element.k = k;
                element.len = len;
                element.q = q;
                element.nnos = nnos; 
                element.face = face;
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
            memb.nnos = 0; 
            memb.face = 0; 
            memb.cord = 0; 
        end
    end
end