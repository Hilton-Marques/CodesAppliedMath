classdef Geometry < handle
    properties
        m_vertices
    end
    methods
        function this = Geometry(vertices)
            this.m_vertices = vertices;
        end
        function coord = inputVertexPosition(this, id)
            coord = this.m_vertices(id,:);
        end
    end
end