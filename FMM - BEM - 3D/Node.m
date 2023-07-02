classdef Node < handle
    properties
        element = [];
        pos = [];
    end
    methods
        function this = Node(Element)
            if nargin > 0
                this.element = Element;
                this.pos = sum(Element.pos)/3;
            end
        end
    end
end