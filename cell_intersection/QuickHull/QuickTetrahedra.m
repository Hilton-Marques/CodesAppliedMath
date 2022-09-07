classdef QuickTetrahedra < handle
    properties
        triangle
        v4
        pts
        halfedges
        faces
        n_true_faces
        triangles = {};
        inc = []
    end
    methods
        function this = QuickTetrahedra(triangle,v4)
            this.triangle = triangle;
            this.v4 = v4;
            this.markVerts();
            this.pts = [triangle.pts;v4.coord];
            this.halfedges = HalfEdge();
            this.halfedges(12) = HalfEdge();
            this.faces = Face();
            this.faces(4) = Face();
            % Create Halfs
            
            % primeiro triangulo ABC
            
            this.halfedges(1).begin_vert = triangle.v1;
            this.halfedges(1).end_vert = triangle.v2;
            this.halfedges(1).face_id = 1;
            this.halfedges(1).next = this.halfedges(2);
            

            this.halfedges(2).begin_vert = triangle.v2;
            this.halfedges(2).end_vert = triangle.v3;
            this.halfedges(2).face_id = 1;
            this.halfedges(2).next = this.halfedges(3);
            
            this.halfedges(3).begin_vert = triangle.v3;
            this.halfedges(3).end_vert = triangle.v1;
            this.halfedges(3).face_id = 1;
            this.halfedges(3).next = this.halfedges(1);
            
            % segundo triangulo BAD            
          
            this.halfedges(4).begin_vert = triangle.v2;
            this.halfedges(4).end_vert = triangle.v1;
            this.halfedges(4).face_id = 2;
            this.halfedges(4).next = this.halfedges(5);
            
            
            this.halfedges(5).begin_vert = triangle.v1;
            this.halfedges(5).end_vert = this.v4;
            this.halfedges(5).face_id = 2;
            this.halfedges(5).next = this.halfedges(6);
            
            this.halfedges(6).begin_vert = this.v4;
            this.halfedges(6).end_vert = triangle.v2;
            this.halfedges(6).face_id = 2;
            this.halfedges(6).next = this.halfedges(4);            
            
            % Terceiro triangulo CBD
            
            this.halfedges(7).begin_vert = triangle.v3;
            this.halfedges(7).end_vert = triangle.v2;          
            this.halfedges(7).face_id = 3;
            this.halfedges(7).next = this.halfedges(8);
            
            this.halfedges(8).begin_vert = triangle.v2;
            this.halfedges(8).end_vert = this.v4;         
            this.halfedges(8).face_id = 3;
            this.halfedges(8).next = this.halfedges(9);
                  
            this.halfedges(9).begin_vert = this.v4;
            this.halfedges(9).end_vert = triangle.v3;         
            this.halfedges(9).face_id = 3;
            this.halfedges(9).next = this.halfedges(7);
            
            % quarto triangulo            
            this.halfedges(10).begin_vert = triangle.v1;
            this.halfedges(10).end_vert = triangle.v3;         
            this.halfedges(10).face_id = 4;
            this.halfedges(10).next = this.halfedges(11);
            
            this.halfedges(11).begin_vert = triangle.v3;
            this.halfedges(11).end_vert = this.v4;         
            this.halfedges(11).face_id = 4;
            this.halfedges(11).next = this.halfedges(12);
            
            this.halfedges(12).begin_vert = this.v4;
            this.halfedges(12).end_vert = triangle.v1;         
            this.halfedges(12).face_id = 4;
            this.halfedges(12).next = this.halfedges(10);    
            
            % faces
            count = 1;
            for i = 1:4
                this.faces(i) = Face(this.halfedges(count));
                count = count + 3;
            end
            %create tiwns
            this.halfedges(1).twin = this.halfedges(4);
            this.halfedges(4).twin = this.halfedges(1);
            this.halfedges(2).twin = this.halfedges(7);
            this.halfedges(7).twin = this.halfedges(2);
            this.halfedges(3).twin = this.halfedges(10);
            this.halfedges(10).twin = this.halfedges(3);
            this.halfedges(5).twin = this.halfedges(12);
            this.halfedges(12).twin = this.halfedges(5);
            this.halfedges(6).twin = this.halfedges(8);
            this.halfedges(8).twin = this.halfedges(6);
            this.halfedges(9).twin = this.halfedges(11);
            this.halfedges(11).twin = this.halfedges(9);
            % add Point
            
            
        end
        function addHed(this)
            this.halfedges(end+1) = HalfEdge();
        end
        function addFace(this,hed)
            this.faces(end+1) = Face(hed);
        end
        function show(this,color) 
            count = 0;
            for f = this.faces
                if ~(f.disabled)
                    this.inc(end+1,1:3) = f.inc;
                    count = count + 1;
                    if nargin == 1
                        color = rand(1,3);
                    end
                    f.show(color);
                end
            end
            this.n_true_faces = count;
        end
        function inc = getIds(this)
            count = 0;
            this.pts = [];
            this.triangles = Face.empty;
            for f = this.faces
                if ~(f.disabled)
                    this.inc(end+1,1:3) = f.inc;
                    this.triangles(end+1) = f;
                    count = count + 1;
                end
            end
            this.n_true_faces = count;
            inc = this.inc;
        end
        function volume = getVolume(this)
            n = size(this.inc,1);
            volume = 0.0;
            for i = 1:n
                f = this.triangles(i);
                qi = f.pts(1,:);
                qj = f.pts(2,:);
                qk = f.pts(3,:);
                volume_i = (dot(qi,cross(qj,qk)));
                volume = volume + volume_i;
            end
            volume = volume/6;
        end
        function markVerts(this)
            this.triangle.v1.Marked = true;
            this.triangle.v2.Marked = true;
            this.triangle.v3.Marked = true;
            this.v4.Marked = true;
        end
    end
end