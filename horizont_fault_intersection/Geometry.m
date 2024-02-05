classdef Geometry < handle
    properties
        m_vertices
    end
    methods
        function this = Geometry(vertices)
            this.m_vertices = vertices;
        end
        function coord = inputVertexPosition(this, id)
            if nargin == 1
                coord = this.m_vertices;
            else
                coord = this.m_vertices(id,:);
            end
        end
    end
end