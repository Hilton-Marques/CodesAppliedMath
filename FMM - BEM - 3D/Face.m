classdef Face < handle
    properties
        q = [];
    end
    methods
        function this = Face(q)
            if nargin > 0
                this.q = q;
            end
        end
    end
end