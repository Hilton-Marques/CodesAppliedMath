classdef Node < handle
    properties
        Element = [];
        pos = [];
    end
    methods
        function this = Node(Element)
            if nargin > 0
                this.Element = Element;
                this.pos = sum(Element.points)/3;
            end
        end
    end
end