classdef Element < handle
    properties
        k = 1;
        len;
        q = [];
        nNos;
        face = [];
        position = [];
        id = 1;
    end
    methods
        function this = Element(nNos,len,face,position,id)
            if nargin > 0
                this.len = len;
                this.nNos = nNos;
                this.position = position;
                this.face = face;
                this.q = face.q;
                this.id = id;
            end
        end
    end
end