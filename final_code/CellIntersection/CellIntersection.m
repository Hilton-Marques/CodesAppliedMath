classdef CellIntersection < handle
    properties
        m_cell_A
        m_cell_B
        m_planes_A
        m_planes_B
        m_convex_A
        m_convex_B
        m_threshold
        m_threshold_volume_abs = 1000
        m_p % interior point
        m_drawer Drawer
        m_gjk_obj GJK
        m_quick_hull QuickHull
    end
    methods
        function [this,volume] = CellIntersection(cell_A, cell_B, threshold)
            this.m_cell_A = cell_A;
            this.m_cell_B = cell_B;
            this.m_threshold = threshold;
            this.m_drawer = Drawer(this);
            volume = this.compute();
        end
        function volume = compute(this)
            volume = 0.0;
            this.m_drawer.showCells();
            this.m_planes_A = this.buildHalfPlanes(this.m_cell_A);
            this.m_planes_B = this.buildHalfPlanes(this.m_cell_B);          
            this.m_convex_A = ConvexCell(this.projCoords(this.m_cell_A,this.m_planes_A));
            this.m_convex_B = ConvexCell(this.projCoords(this.m_cell_B,this.m_planes_B));
            this.m_convex_A.show('red');
            this.m_convex_B.show('blue');
            this.m_drawer.showConvexApproximation();      
            [isIntersecting, p, planes] = this.isIntersecting(this.m_cell_A , this.m_cell_B);
            volume_i = this.getIntersectionVolume(planes, p);
            [err_concave, exact,convex_volume] = this.errConcave(this.m_cell_A );
            %isIntersecting = true;

            %this.m_drawer.showGJK();
            %if isIntersecting
            if true
                %volume = this.getIntersectionVolume([this.m_planes_A;this.m_planes_B], this.m_p);
                new_cells_A = {};
                new_cells_A = this.recSubdivide(this.m_cell_A,new_cells_A);
                new_cells_B = {};
                new_cells_B = this.recSubdivide(this.m_cell_B,new_cells_B);   
                for j = 1:size(new_cells_B,2)
                    cell_B = new_cells_B{j};
                    %this.m_drawer.showSurface(cell_B, this.m_drawer.m_blue);
                end
                for i = 1:size(new_cells_A,2)
                    cell_A = new_cells_A{i};
                    this.m_drawer.showSurface(cell_A, this.m_drawer.m_red);
                    for j = 1:size(new_cells_B,2)
                        cell_B = new_cells_B{j};
                        %this.m_drawer.showParentChildren(cell_A , {cell_B});
                        [isIntersecting, p, planes] = this.isIntersecting(cell_A , cell_B);
                        if isIntersecting
                            this.m_drawer.showParentChildren(cell_A , {cell_B});
                            %this.m_drawer.showSurface(cell_B, this.m_drawer.m_blue);
                            %this.showChildGJK(p,cell_A,cell_B);
                            volume_i = this.getIntersectionVolume(planes, p);
                            volume = volume + volume_i;
                        end
                    end
                end
                 U(1,:) = reshape(this.m_cell_A',24,1)';
                 U(2,:) = reshape(this.m_cell_B',24,1)';
                 this.m_drawer.Shadow(U,1);
                 this.m_drawer.setVolumeTitle(volume);
            end
        end
        function showChildGJK(this,p,cell_A,cell_B)
            axis tight
            plot3(p(1),p(2),p(3), 'o','MarkerSize',11,'MarkerFaceColor','black','Color','black')
            U(1,:) = reshape(cell_A',24,1)';
            U(2,:) = reshape(cell_B',24,1)';
            this.m_drawer.Shadow(U,1);
            exportgraphics(gca,strcat('gjk','.pdf'),'Resolution',2000)
        end
        function [bool, p, planes] = isIntersecting(this, cell_A,cell_B)
            planes_A = this.buildHalfPlanes(cell_A);
            planes_B = this.buildHalfPlanes(cell_B);
            convex_A = ConvexCell(this.projCoords(cell_A,planes_A,false));
            convex_B = ConvexCell(this.projCoords(cell_B,planes_B,false));
            planes = [planes_A;planes_B];
            this.m_gjk_obj = GJK(convex_A,convex_B);
            [bool, p] = this.m_gjk_obj.compute(planes);
        end
        function res = getCoords(this)
            res = [this.m_cell_A;this.m_cell_B];
        end
        function new_cells = recSubdivide(this, pts, new_cells)
            [err_concave, total_volume] = this.errConcave(pts);
            if (err_concave < this.m_threshold || total_volume < this.m_threshold_volume_abs)  
                new_cells{end+1} = pts;
            else
                children_cells = this.childrenCells(pts);
                for i = 1:4
                    new_cells = this.recSubdivide(children_cells{i}, new_cells);
                end
            end
        end
        function volume = getIntersectionVolume(this, planes, center)
            coord_dual = [];
            A = planes(:,1:3);
            b = -planes(:,4);
            for i = 1:size(planes,1)
                n = A(i,1:3);
                d = b(i);
                new_d = d - dot(center,n)
                if (new_d == 0)
                    continue;
                end
                coord_dual(end+1,1:3) = ( 1/new_d ) * n;
            end
            %[k1,vol] = convhull(coord_dual);
            [~,k1] = QuickHull(coord_dual);
            dualPoints = zeros(size(k1,1),3);
            count = 1;
            for i = 1:size(k1,1)
                inc = k1(i,:);
                p1 = coord_dual(inc(1),:);
                p2 = coord_dual(inc(2),:);
                p3 = coord_dual(inc(3),:);
                normal = cross(p2 - p1, p3 - p1);
                normal = normal / norm(normal);
                d = dot(normal,p1);
                if (d == 0)
                    continue;
                end
                dualPoints(count,:) = (1/d)* normal + center  ;
                count = count + 1;
            end
            [k2,volume] = convhull(dualPoints);
            %[~,k2,volume] = QuickHull(dualPoints);
            trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','magenta','EdgeAlpha',1.0);
        end
    end
    methods (Static)
        %geometry utils
       function planes = buildHalfPlanes(pts)
            faces = cell(6,1);
            faces(1) = {[ pts(1,:) ; pts(2,:) ; pts(4,:) ; pts(3,:)] };
            faces(2) = {[ pts(5,:) ; pts(7,:) ; pts(8,:) ; pts(6,:)] };
            faces(3) = {[ pts(2,:) ; pts(6,:) ; pts(8,:) ; pts(4,:)] };
            faces(4) = {[ pts(1,:) ; pts(3,:) ; pts(7,:) ; pts(5,:)] };
            faces(5) = {[ pts(1,:) ; pts(5,:) ; pts(6,:) ; pts(2,:)] };
            faces(6) = {[ pts(3,:) ; pts(4,:) ; pts(7,:) ; pts(7,:)] };
            planes = zeros(6,4);
            for i = 1:6
                face_pts = faces{i};
                planes(i,:) = CellIntersection.fittingPlane(face_pts);
            end
            opp_ids = [2,1,4,3,6,5];
            for i = 1:6
                plane = planes(i,:);
                normal = plane(1:3);
                %correct plane
                if (norm(normal) == 0)
                    opp_plane = planes(opp_ids(i),:);
                    normal = opp_plane(1:3);
                    if (norm(normal) == 0)
                        % não tem o que fazer
                    else
                        max_dist = - realmax;
                        id_far = -1;
                        plane_pts = faces{i};
                        pn = - opp_plane(1:3) * opp_plane(4);
                        for j = 1:4
                            pi = plane_pts(j,:) - pn;
                            dist = abs(dot(pi, normal));
                            if dist > max_dist
                                id_far = j;
                                max_dist = dist;
                            end
                        end
                        planes(i,:) = [-normal, dot(plane_pts(id_far,:),normal)];
                    end
                end
            end
        end
        function planes = buildHalfPlanesTopBottom(pts)
            faces = cell(6,1);
            faces(1) = {[ pts(1,:) ; pts(2,:) ; pts(4,:) ; pts(3,:)] };
            faces(2) = {[ pts(5,:) ; pts(7,:) ; pts(8,:) ; pts(6,:)] };
            faces(3) = {[ pts(2,:) ; pts(6,:) ; pts(8,:) ; pts(4,:)] };
            faces(4) = {[ pts(1,:) ; pts(3,:) ; pts(7,:) ; pts(5,:)] };
            faces(5) = {[ pts(1,:) ; pts(5,:) ; pts(6,:) ; pts(2,:)] };
            faces(6) = {[ pts(3,:) ; pts(4,:) ; pts(7,:) ; pts(7,:)] };
            %planes = zeros(6,4);
            planes = zeros(8,4);
            %planes(1:2,:) = CellIntersection.fittingTopBottom(faces{1},1);
            pts = faces{1};
            pt = (pts(2,:) + pts(4,:)) * 0.5;
            a = CellIntersection.getNormal(pts([1,2,4],:));
            b = CellIntersection.getNormal(pts([2,3,4],:));
            planes(1:2,:) = [[a,-dot(a, pt)];[b,-dot(b,pt)]];
            pts = faces{2};
            pt = (pts(1,:) + pts(3,:)) * 0.5;
            a = CellIntersection.getNormal(pts([1,2,3],:));
            b = CellIntersection.getNormal(pts([1,3,4],:));
            planes(3:4,:) = [[a,-dot(a, pt)];[b,-dot(b,pt)]];
            %for i = 3:6
            for i = 3:6
                face_pts = faces{i};
                planes(i+2,:) = CellIntersection.fittingPlane(face_pts);
            end
            opp_ids = [2,1,4,3,6,5];
            for i = 3:6
                %plane = planes(i,:);
                plane = planes(i+2,:);
                normal = plane(1:3);
                %correct plane
                if (norm(normal) == 0)
                    opp_plane = planes(opp_ids(i)+2,:);
                    normal = opp_plane(1:3);
                    if (norm(normal) == 0)
                        % não tem o que fazer
                    else
                        max_dist = - realmax;
                        id_far = -1;
                        plane_pts = faces{i};
                        pn = - opp_plane(1:3) * opp_plane(4);
                        for j = 1:4
                            pi = plane_pts(j,:) - pn;
                            dist = abs(dot(pi, normal));
                            if dist > max_dist
                                id_far = j;
                                max_dist = dist;
                            end
                        end
                        planes(i+2,:) = [-normal, dot(plane_pts(id_far,:),normal)];
                    end
                end
            end
        end
        function [o,d] = intersect3D(n,m)
            z = cross(m,n);
            m_perp = cross(z,m);
            c = dot(n,n);
            d = dot(m,n);
            e = dot(m_perp,n);
            lam2 = (c-d)/e;
            o = m + lam2*m_perp;
            lam = 10000;
            p1 = o + lam*z/norm(z);
            p2 = o - lam*z/norm(z);
            d = z/norm(z);
            %line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'linewidth',3,'color','black');
            check1 = dot(o-n,n);
            check2 = dot(o-m,m);
        end
        function p = rayPlaneIntersection(o,ray,n)
            d = n(4);
            n_vec = n(1:3);
            d_i = dot(n_vec,o) + d;
            lam = - d_i / dot(n_vec,ray);
            p = o + lam*ray;
        end
        function planes = fittingTopBottom(pts,sign)
            pt = (pts(2,:) + pts(4,:)) * 0.5;
            a = sign*CellIntersection.getNormal(pts([1,2,4],:));
            b = sign*CellIntersection.getNormal(pts([2,3,4],:));
            planes = [[a,-dot(a, pt)];[b,-dot(b,pt)]];
        end
        function n = getNormal(pts)
            p0 = pts(1,:);
            p1 = pts(2,:);
            p2 = pts(3,:);            
            n = cross(p1 - p0,p2 - p0);
            len = norm(n);
            if len == 0
                n = zeros(3,1);
            else
                n = n/len;
            end
        end
        function plane = fittingPlane(pts)

            centroide = mean(pts);
            normal = zeros(1,3);
            plane = zeros(1,4);
            %k = convhull(pts(:,1),pts(:,2),pts(:,3));

            %a = getNormal(pts([1,2,3],:));
            % b = getNormal(pts([1,3,4],:));
            % if dot(a,b) < 0
            %     b = -b;
            % end
            % normal = a + b;
            % for i = 1:size(k,2)
            %     normal = normal + getNormal(pts(k(i,:),:));
            % end
            for i = 1:4
                p0 = pts(i,:);
                p1 = pts(mod(i,4) + 1 , :);
                p2 = pts(mod(i+1,4) + 1 , :);
                normal_i = cross(p1-p0,p2-p0);
                len = norm(normal_i);
                if len == 0
                    continue
                end
                normal_i = normal_i / len;
                normal = normal + normal_i;
            end
            if norm(normal) == 0
                return
            end
            normal = normal/norm(normal);
            plane = [normal,-dot(normal, centroide)];
        end
        function new_coords = projCoords(pts,planes,flag)
%             new_coords = pts;
%             return
            if nargin == 2
                flag = false;
            end
            new_coords = zeros(8,3);

            ids_pts = [1,2,4,3,5,6,8,7];
            pts = pts(ids_pts,:);   % put on circular order

            ids_planes = [4,5,3,6,1,2]; % lef front right back top bottom
            planes = planes(ids_planes,:);

            [o,d] = CellIntersection.intersect3D(-planes(6,4)*planes(6,1:3),-planes(5,4)*planes(5,1:3));
            edges_len = zeros(1,4);
            count = 0;
            if flag
                lam = 10000;
                p1 = o + lam*d;
                p2 = o - lam*d;
                line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'linewidth',3,'color','black');
            end
            for i = 1:4
                p_i = pts(i,:);
                p_j = pts(i+4,:);
                edge = p_i - p_j;
                edges_len(i) = dot(edge,edge);
                if (edges_len(i) == 0)
                    count = count + 1;
                    edge_dir = cross(planes(i,1:3),planes(mod(i,4)+1,1:3));
                    p_i = p_j + edge_dir;
                    p_j = p_j - edge_dir;
                    edge = p_i - p_j;
                end
                p_i = CellIntersection.rayPlaneIntersection(p_i,edge,planes(5,:));
                p_j = CellIntersection.rayPlaneIntersection(p_j,edge,planes(6,:));
                new_edge = p_i - p_j;
                if dot(edge,new_edge) < 0
                    p_j = CellIntersection.rayPlaneIntersection(o,d,planes(i,:));
                    p_i = CellIntersection.rayPlaneIntersection(o,d,planes(mod(i,4)+1,:));
                end
                edge = p_i - p_j;
                len = norm(edge);
                edge = edge/len;
                new_coords(i,:) = p_i - edge*1e-2;
                new_coords(i + 4,:) = p_j + edge*1e-2;
            end
            %check degenerateFaces
            if count == 2
                for i = 1:4
                    ei = edges_len(i);
                    j = mod(i,4)+1;
                    ej = edges_len(j);
                    if (ei == 0 ) && (ej == 0)
                        %denerate face
                        p_i = ReservoirIntersection.rayPlaneIntersection(o,d,planes(i,:));
                        p_j = ReservoirIntersection.rayPlaneIntersection(o,d,planes(mod(j,4)+1,:));
                        new_coords(i,:) = p_i;
                        new_coords(i + 4,:) = p_i;
                        new_coords(j,:) = p_j;
                        new_coords(j + 4,:) = p_j;
                    end
                end
            end
        end
        function children_cells = childrenCells(pts)
            top = [pts(1,:);pts(2,:);pts(4,:); pts(3,:)];
            bottom = [pts(5,:);pts(6,:);pts(8,:); pts(7,:)];
            normals = CellIntersection.createNormals(top);

            new_pts_top = CellIntersection.subdivide(top);
            new_pts_bottom = CellIntersection.subdivide(bottom);
            solid1 = [[top(1,:) ; new_pts_top(1,:); new_pts_top(5,:);new_pts_top(4,:)];...
                [bottom(1,:) ; new_pts_bottom(1,:); new_pts_bottom(5,:);new_pts_bottom(4,:)]];


            solid2 = [[new_pts_top(1,:); top(2,:); new_pts_top(2,:); new_pts_top(5,:)];...
                [new_pts_bottom(1,:); bottom(2,:); new_pts_bottom(2,:); new_pts_bottom(5,:)]];


            solid3 = [[new_pts_top(5,:); new_pts_top(2,:); top(3,:);new_pts_top(3,:)];...
                [new_pts_bottom(5,:); new_pts_bottom(2,:);bottom(3,:); new_pts_bottom(3,:)]];

            solid4 = [[new_pts_top(4,:); new_pts_top(5,:);new_pts_top(3,:);top(4,:)];...
                [new_pts_bottom(4,:); new_pts_bottom(5,:);new_pts_bottom(3,:);bottom(4,:)]];

            ids = [1,2,4,3,5,6,8,7];
            solid1 = solid1(ids,:);
            solid2 = solid2(ids,:);
            solid3 = solid3(ids,:);
            solid4 = solid4(ids,:);
            children_cells = {solid1, solid2, solid3, solid4};
        end
        function normals = createNormals(face)
            A = face(1,:)';
            B = face(2,:)';
            C = face(3,:)';
            D = face(4,:)';
            L1 = B - A;
            L2 = C - B;
            L3 = D - C;
            L4 = A - D;
            normals = zeros(4,3);
            normals(1,:) = -L1 + L4;
            normals(2,:) = L1 - L2;
            normals(3,:) = -L3 + L2;
            normals(4,:) = L3 - L4;
            normals(1,:) = normals(1,:)/norm(normals(1,:));
            normals(2,:) = normals(2,:)/norm(normals(2,:));
            normals(3,:) = normals(3,:)/norm(normals(3,:));
            normals(4,:) = normals(4,:)/norm(normals(4,:));
        end
        function [ratio,volume_total,volume_max_convex] = errConcave(pts)
            planes{1} = {pts(1,:),pts(2,:),pts(4,:), pts(3,:)};
            planes{2} = {pts(5,:),pts(7,:),pts(8,:), pts(6,:)};
            planes{3} = {pts(1,:),pts(3,:),pts(7,:), pts(5,:)};
            planes{4} = {pts(2,:),pts(6,:),pts(8,:), pts(4,:)};
            planes{5} = {pts(1,:),pts(5,:),pts(6,:), pts(2,:)};
            planes{6} = {pts(3,:),pts(4,:),pts(8,:), pts(7,:)};
            R = sum(pts,1)/8;
            R = R';
            volume_total = 0.0;
            volume_concave = 0.0;
            volume_convex = 0.0;
            volume_max_convex = 0.0;
            for i = 1:6
                coords = planes{i};
                A = coords{1}';
                B = coords{2}';
                C = coords{3}';
                D = coords{4}';
                volume = (1/6)*dot(A-R,cross((B-D),(C-A))) + (1/12)*dot(C-A,cross((D-A),(B-A)));
                conv = (1/6)*dot(A-R,cross((B-D),(C-A)));
                conc = (1/12)*dot(C-A,cross((D-A),(B-A)));
                if  conc > 0
                    volume_max_convex =  volume_max_convex + (1/6)*dot(B-R,cross((B-D),(C-A)));
                else
                    volume_max_convex =  volume_max_convex + conv;
                end
                
                volume_total = volume_total + volume;
                volume_concave = volume_concave + abs(((1/12)*dot(C-A,cross((D-A),(B-A)))));
                %volume4 = (1/6)*dot(B-R,cross((B-D),(C-A))) - (1/12)*dot(B-A,cross((C-A), (D-A))) + (1/6)*dot(B-A,cross((C-A),(B-D)));
                volume7 = (1/6)*dot(B-R,cross((B-D),(C-A))) + (1/12)*dot(D-B,cross((A-B), (C-D))) ;
                tetra = [R';A';B';C';D'];
                %volume2 = getVol(tetra,[2,3,4;1,3,2;1,4,3;1,2,4]) + getVol(tetra,[2,4,5;1,4,2;1,5,4;1,2,5]);
                %volume3 = getVol(tetra,[2,3,5;1,3,2;1,5,3;1,2,5]) + getVol(tetra,[3,4,5;1,4,3;1,5,4;1,3,5]);
                volume_convex = volume_convex + (1/6)*dot(A-R,cross((B-D),(C-A)));
            end
            %ratio = (volume_total - volume_convex) / volume_total;
            ratio = volume_concave / volume_total;
        end
        function new_pts = subdivide(face)
            A = face(1,:)';
            B = face(2,:)';
            C = face(3,:)';
            D = face(4,:)';
            new_pts = zeros(5,3);
            ksi_midle = [0.5,1,0.5,0,0.5];
            eta_midle = [0,0.5,1.0,0.5,0.5];
            x = A(1) + ksi_midle.*(B(1) - A(1)) + eta_midle.*(D(1) - A(1)) + ksi_midle.*eta_midle.*( (C(1)-A(1)) - (B(1)-A(1)) - (D(1) - A(1)));
            y = A(2) + ksi_midle.*(B(2) - A(2)) + eta_midle.*(D(2) - A(2)) + ksi_midle.*eta_midle.*((C(2)-A(2)) - (B(2)-A(2)) - (D(2) - A(2)));
            z = A(3) + ksi_midle.*(B(3) - A(3)) + eta_midle.*(D(3) - A(3)) + ksi_midle.*eta_midle.*((C(3)-A(3)) - (B(3)-A(3)) - (D(3) - A(3)));
            for i = 1:5
                new_pts(i,:) = [x(i),y(i),z(i)];
            end
        end
    end
end