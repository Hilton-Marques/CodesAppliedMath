classdef QuickHull < handle
    properties
        verts;
        n;
        faceData = struct('face', Face.empty,'hed',HalfEdge.empty);
    end
    methods
        function [this,k,volume] = QuickHull(pts)
            n = size(pts,1);
            verts(n) = Verts();
            for i = 1:n
                p = pts(i,:);
                verts(i) = Verts(i,p);
            end
            this.verts = verts;
            this.n = length(verts);
            tetrahedra = this.getPoly();
            k = tetrahedra.getIds();            
            volume = tetrahedra.getVolume();
        end

        function tetrahedra = getPoly(this)
            % find extremal points
            extreme_values = [this.verts(1).coord(1),this.verts(1).coord(1),...
                this.verts(1).coord(2),this.verts(1).coord(2),...
                this.verts(1).coord(3),this.verts(1).coord(3)];
            extremal_indices = ones(1,6);
            for i = 1:this.n
                vi = this.verts(i);
                if vi.coord(1) < extreme_values(1)
                    extreme_values(1) = vi.coord(1);
                    extremal_indices(1) = vi.id;
                elseif vi.coord(1) > extreme_values(2)
                    extreme_values(2) = vi.coord(1);
                    extremal_indices(2) = vi.id;
                end
                if vi.coord(2) < extreme_values(3)
                    extreme_values(3) = vi.coord(2);
                    extremal_indices(3) = vi.id;
                elseif vi.coord(2) > extreme_values(4)
                    extreme_values(4) = vi.coord(2);
                    extremal_indices(4) = vi.id;
                end
                if vi.coord(3) < extreme_values(5)
                    extreme_values(5) = vi.coord(3);
                    extremal_indices(5) = vi.id;
                elseif vi.coord(3) > extreme_values(6)
                    extreme_values(6) = vi.coord(3);
                    extremal_indices(6) = vi.id;
                end
            end
            %find extremal edge
            max_dist = -1;
            for i = 1:6
               pi = this.verts(extremal_indices(i)).coord;
                for j = i:6
                    pj = this.verts(extremal_indices(j)).coord;
                    u = pi - pj;
                    d_i = dot(u,u);
                    if (d_i > max_dist)
                        max_dist = d_i;
                        edge = [extremal_indices(i),extremal_indices(j)];
                    end
                end
            end
            long_ray = Ray(this.verts(edge(1)),this.verts(edge(2)));
            %find extremal triangle
            max_dist = -1;
            selected_id = -1;
            dists = [];
            for i = 1:this.n
                vi = this.verts(i);
                dist_i = long_ray.getSquaredDistanceBetweenPoint(vi.coord);
                dists(end+1) = dist_i;
                if dist_i > max_dist
                    max_dist = dist_i;
                    selected_id = i;
                end                
            end
            base_triangle = Triangle(this.verts(edge(1)), this.verts(edge(2)), this.verts(selected_id));
            %find extremal tetrahedra
            max_dist = -1;
            selected_id = -1;
            for i = 1:this.n
                vi = this.verts(i);
                dist_signed = base_triangle.getDistanceFromPoint(vi.coord);
                dist_i = abs(dist_signed);
                if dist_i > max_dist
                    max_dist = dist_i;
                    selected_id = i;
                    signed = dist_signed;
                end
            end
            % se signed é positivo indica que a normal está dentro do
            % tetraedro
            if signed > 0
                base_triangle = Triangle(base_triangle.v1,base_triangle.v3, base_triangle.v2);
            end
            tetrahedra = QuickTetrahedra(base_triangle,this.verts(selected_id));
            this.buildEyeTetrahedra(tetrahedra);
            %tetrahedra.show();
            %% inicializa a pilha de faces que possuem pontos a vista
            faceList = Face.empty;
            for f = tetrahedra.faces
                if ~isempty(f.pointsOnPositiveSide)
                    faceList(end + 1) = f;
                end
            end
            % Loop princial: processe as faces até a lista esteja vazia
            iter = 0;
            while ~isempty(faceList)
                iter = iter + 1;
                face_iter = faceList(end);
                faceList(end) = [];
                % get most distant point
                far_point_iter = face_iter.most_distant_point;
                %Construa o horizonte
                horizonHeds = HalfEdge.empty;
                possiblyVisibleFaces = this.faceData();
                visibleFaces = Face.empty;
                %first tri
                possiblyVisibleFaces.face = face_iter;possiblyVisibleFaces.hed = HalfEdge();                
                % collect o horizonte para este ponto em especifico
                while ~isempty(possiblyVisibleFaces)
                    possible_visible_face = possiblyVisibleFaces(end);
                    face = possible_visible_face.face;
                    entered_hed = possible_visible_face.hed;
                    possiblyVisibleFaces(end) = [];                    
                    if face.was_visited_iter == iter
                        if (face.isVisibleFaceOnCurrentIter == true)
                            continue;
                        end
                    else
                        face.was_visited_iter = iter;
                        % é visivel?
                        a = face.getDistanceFromPoint(far_point_iter.coord);
                        if (face.getDistanceFromPoint(far_point_iter.coord) > 0)
                            face.isVisibleFaceOnCurrentIter = true;
                            visibleFaces(end+1) = face;                            
                            for hed = face.heds
                                %busque os vizinhos, porém evite o
                                %triângulo inicial 
                                if (hed.twin ~= entered_hed)
                                    new_possible_face = this.faceData();
                                    new_possible_face.face = hed.twin.face;
                                    new_possible_face.hed = hed;
                                    possiblyVisibleFaces(end+1) = new_possible_face;
                                end
                            end
                            continue;
                        end
                    end
                    % se a face foi visitada e não é visivel então a aresta
                    % de contao faz parte do horizonte
                    %entered_hed.show();
                    horizonHeds(end+1) = entered_hed;
                end
                %this.plotHorizon(horizonHeds,face_iter.color);
                % all visible faces should be disabled
                for face = visibleFaces
                    face.disabled = true;
                end
                n_horizons = length(horizonHeds);
                %reorder horizons in circular order
                for i = 1:n_horizons
                    hed_i_id = horizonHeds(i).end_vert.id;
                    for j = i+1:n_horizons
                        hed_j_id = horizonHeds(j).begin_vert.id;
                        if hed_j_id == hed_i_id
                            temp = horizonHeds(j);
                            horizonHeds(j) = horizonHeds(i+1);
                            horizonHeds(i+1) = temp;
                        end
                    end
                end
                %Update previously the number of new heds
                n_horizons = length(horizonHeds);
                new_heds_id = zeros(1,n_horizons);
                new_faces_id = zeros(1,n_horizons);
                for i = 1:2*n_horizons
                    new_heds_id(i) = length(tetrahedra.halfedges) + 1;
                    tetrahedra.addHed();
                end
                for i = 1:n_horizons
                    AB = horizonHeds(i);
                    A = AB.begin_vert;
                    B = AB.end_vert;
                    C = far_point_iter;                    
                    new_face_id = length(tetrahedra.faces) + 1;
                    new_faces_id(i) = new_face_id;
                    %Create BC 
                    BC = tetrahedra.halfedges(new_heds_id(2*i));                    
                    BC.begin_vert = B;
                    BC.end_vert = C;
                    BC.face_id = new_face_id;
                    %Create CA
                    CA = tetrahedra.halfedges(new_heds_id(2*i - 1));                    
                    CA.begin_vert = C;
                    CA.end_vert = A;
                    CA.face_id = new_face_id;
                    %Update next
                    AB.next = BC;
                    BC.next = CA;
                    CA.next = AB;
                    %Create new Face
                    tetrahedra.addFace(AB);
                    % update oposite
                    BC.twin = tetrahedra.halfedges(new_heds_id(mod((i)*2 , 2*n_horizons) + 1)); 
                    CA.twin = tetrahedra.halfedges(new_heds_id(mod( 2*(i-1) - 1 , 2*n_horizons) + 1)); 
                end
                %update cloud set
                for face = visibleFaces
                    for pt = face.pointsOnPositiveSide
                        if pt.id == far_point_iter.id
                            continue
                        end
                        % paracada nova face
                        for i = 1:n_horizons
                            new_face = tetrahedra.faces(new_faces_id(i));
                            if (new_face.addPoint(pt))
                                break;
                            end                                
                        end
                    end
                end
                %update face stack
                for i = 1:n_horizons
                    face = tetrahedra.faces(new_faces_id(i));
                    if ~isempty(face.pointsOnPositiveSide)
                        faceList(end+1) = face;
                    end
                end
            end
        end

        function buildEyeTetrahedra(this,tetrahedra)
            for v = this.verts
                if ~v.Marked
                    for face = tetrahedra.faces
                        if face.addPoint(v)
                            break;
                        end
                    end
                end
            end
        end
    end
    methods (Static)        
        function area2 = squaredArea(p0,p1,p2)
            u = p1 - p0;
            v = p2 - p0;
            area = cross(u,v);
            area2 = dot(area,area);
        end
        function volume = signedVolume(p0,p1,p2,p3)
            u = p1 - p0;
            v = p2 - p0;
            z = p3 - p0;
            volume = det([u;v;z]);
            if abs(volume) < 0.0001
                volume = 0;
            end
        end
        function plotHorizon(horizons,color)
            for hed = horizons
                hed.show(color);
            end
        end
    end

    
end