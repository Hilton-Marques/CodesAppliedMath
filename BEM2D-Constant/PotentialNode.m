%
%% Class definition
classdef PotentialNode < handle
    %%
    % <https://www.mathworks.com/help/matlab/ref/handle-class.html
    % See documentation on *handle* super-class>.
    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        cord = []; %cord per nodes
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function node = CrossNode(cord)
            if (nargin > 0)
                 node.cord = cord;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of a CrossNode object.
        function node = clean(node)
                node.cord = 0;
        end
    end
end