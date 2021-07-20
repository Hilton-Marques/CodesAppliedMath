classdef Node < handle
    properties
        position = [];
        element = [];
    end
    methods
        function this = Node(pos,element)
            if nargin > 0
                this.position = pos;
                this.element = element;
            end
        end
    end
end