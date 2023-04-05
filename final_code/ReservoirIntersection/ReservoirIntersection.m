classdef ReservoirIntersection < handle
    properties
        m_model_A
        m_model_B
        m_drawer
        m_threshold
        m_bb
        m_axis =  1.0e+03*[2.5766,3.1766,-1.1760,-0.5408,-3.2680,-3.1949];
        m_ids
    end
    methods
        function this = ReservoirIntersection(model_A,model_B,threshold,ids_A,ids_B)
            this.m_model_A = model_A;
            this.m_model_B = model_B;
            %this.m_drawer = Drawer(this);
            this.m_bb = BBIntersection(model_A,ids_A);
            this.m_threshold = threshold;  
            this.m_ids = ids_A;
            %err = this.CheckSize();
            volume = this.Intersection();
        end
        function err = CheckSize(this)
            e = zeros(size(this.m_model_A,1),1);
            volume_tot = 0.0;
            for i = 1:size(this.m_model_A,1)                 
                cell_A = reshape(this.m_model_A(i,:),3,8)';                
                new_cells_A = {};
                new_cells_A = this.recSubdivide(cell_A,new_cells_A);
                volume = 0.0;
                for child_A = new_cells_A
                    %this.m_drawer.showSurface(child_A{1}, this.m_drawer.m_blue)
                    [isIntersecting, p, planes] = this.isIntersecting(child_A{1}, child_A{1});
                    %plot3(p(1),p(2),p(3),'o','markersize',10);
                    if isIntersecting
                        volume_i = this.getIntersectionVolume(planes, p);
                        volume = volume + volume_i;
                    end
                end
                [ratio,exact] = this.errConcave(cell_A);
                e(i) = abs(exact - volume);
                volume_tot = volume_tot + exact;
            end
            f = sort(e,'descend');
            err = sum(e)/volume_tot;
        end
        function volume = Intersection(this)
            volume = 0;
            %for i = 1:size(model_A,2)
            for i = 30475:size(this.m_model_A,1)                
                cell_A = reshape(this.m_model_A(i,:),3,8)';
                [cells,bins_ids] = this.m_bb.FindGroup(this.m_model_A(i,:));
                new_cells_A = {};
                new_cells_A = this.recSubdivide(cell_A,new_cells_A);
                for j = 1:length(cells)
                    cell_B = reshape(this.m_model_B(cells(j),:),3,8)';
                    new_cells_B = {};
                    new_cells_B = this.recSubdivide(cell_B,new_cells_B);
                    for child_A = new_cells_A
                        %this.m_drawer.showSurface(child_A{1}, this.m_drawer.m_red)
                        for child_B = new_cells_B
                            %this.m_drawer.showSurface(child_B{1}, this.m_drawer.m_blue)
                            [isIntersecting, p, planes] = this.isIntersecting(child_A{1}, child_B{1});
                            if isIntersecting
                                volume_i = this.getIntersectionVolume(planes, p);
                                volume = volume + volume_i;
                            end
                        end
                    end
                end
            end
        end
        function [bool, p, planes] = isIntersecting(this, cell_A,cell_B)
            planes_A = this.buildHalfPlanes(cell_A);
            planes_B = this.buildHalfPlanes(cell_B);
            convex_A = ConvexCell(this.projCoords(cell_A,planes_A,false));
            convex_B = ConvexCell(this.projCoords(cell_B,planes_B,false));
            planes = [planes_A;planes_B];
            gjk_obj = GJK(convex_A,convex_B);
            [bool, p] = gjk_obj.compute(planes);
        end
        function volume = getIntersectionVolume(this, planes, center)
            coord_dual = [];
            A = planes(:,1:3);
            b = -planes(:,4);
            for i = 1:size(planes,1)
                n = A(i,1:3);
                d = b(i);
                new_d = d - dot(center,n);
                if (new_d == 0)
                    continue;
                end
                coord_dual(end+1,1:3) = ( 1/new_d ) * n;
            end
            [k1,vol] = convhull(coord_dual);
            %[~,k1] = QuickHull(coord_dual);
            dualPoints = zeros(size(k1,1),3);
            count = 1;
            for i = 1:size(k1,1)
                inc = k1(i,:);
                p1 = coord_dual(inc(1),:);
                p2 = coord_dual(inc(2),:);
                p3 = coord_dual(inc(3),:);
                normal = cross(p2 - p1, p3 - p1);
                area = norm(normal);

                normal = normal / norm(normal);
                d = dot(normal,p1);
                if (abs(d) < 1e-8)
                    continue;
                end
                dualPoints(count,:) = (1/d)* normal + center  ;
                count = count + 1;
            end
            [k2,volume] = convhull(dualPoints(1:count-1,1:3));
            %[~,k2,volume] = QuickHull(dualPoints);
            %trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','magenta','EdgeAlpha',0.0);
        end
        function new_cells = recSubdivide(this, pts, new_cells)
            [err_concave, total_volume] = this.errConcave(pts);
            if (err_concave < this.m_threshold || total_volume < 1000)
                new_cells{end+1} = pts;
            else
                children_cells = this.childrenCells(pts);
                for i = 1:4
                    new_cells = this.recSubdivide(children_cells{i}, new_cells);
                end
            end
        end
        function drawNeigh(this,i,bins_ids,cells)
            this.m_drawer.showCellss(this.m_model_A(i,:),this.m_drawer.m_blue,1.0);
            this.m_drawer.showGrid(bins_ids,0.8);
            this.m_drawer.showCellss(this.m_model_B(cells,:),this.m_drawer.m_red,0.6);
            axis vis3d;
            camproj('persp');
            T = get(gca,'tightinset');
            set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
            this.m_drawer.Shadow(this.m_model_B(cells,:));
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
            faces(6) = {[ pts(3,:) ; pts(4,:) ; pts(8,:) ; pts(7,:)] };
            planes = zeros(6,4);
            for i = 1:6
                face_pts = faces{i};
                planes(i,:) = ReservoirIntersection.fittingPlane(face_pts);
            end
            opp_ids = [2,1,4,3,6,5];
            del_ids = [];
            for i = 1:6
                plane = planes(i,:);
                normal = plane(1:3);
                %correct plane
                if (norm(normal) == 0)
                    opp_plane = planes(opp_ids(i),:);
                    normal = opp_plane(1:3);
                    if (norm(normal) == 0)
                        del_ids(end+1) = opp_ids(i);
                        del_ids(end+1) = i;
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
        function plane = fittingPlane(pts)
            centroide = mean(pts);
            normal = zeros(1,3);
            plane = zeros(1,4);
            %k = convhull(pts(:,1),pts(:,2),pts(:,3))
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
                if len < 1e-8
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
            if nargin == 2
                flag = true;
            end
            new_coords = zeros(8,3);
            pts_old = pts;
            ids_pts = [1,2,4,3,5,6,8,7];
            pts = pts(ids_pts,:);   % put on circular order
            pts_old_order = pts;
            ids_planes = [4,5,3,6,1,2]; % lef front right back top bottom
            planes = planes(ids_planes,:);
            inverted = -1;
            [o,d] = ReservoirIntersection.intersect3D(-planes(6,4)*planes(6,1:3),-planes(5,4)*planes(5,1:3));
            edges_len = zeros(1,4);
            if flag
                lam = 10000;
                p1 = o + lam*d;
                p2 = o - lam*d;
                line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'linewidth',3,'color','black');
            end
            count = 0;
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
                p_i = ReservoirIntersection.rayPlaneIntersection(p_i,edge,planes(5,:));
                p_j = ReservoirIntersection.rayPlaneIntersection(p_j,edge,planes(6,:));
                new_edge = p_i - p_j;
                if dot(edge,new_edge) < 0
                    p_j = ReservoirIntersection.rayPlaneIntersection(o,d,planes(i,:));
                    p_i = ReservoirIntersection.rayPlaneIntersection(o,d,planes(mod(i,4)+1,:));
                    inverted = i;
                end
                edge = p_i - p_j;
                len = norm(edge);
                edge = edge/len;
                new_coords(i,:) = p_i;
                new_coords(i + 4,:) = p_j;
                %plot3(p_i(1),p_i(2),p_i(3),'o','MarkerSize',10);
                %plot3(p_j(1),p_j(2),p_j(3),'o','MarkerSize',10);
            end
            %check degenerateFaces
            if count == 2  
                for i = 1:4
                    ei = edges_len(i);
                    j = mod(i,4)+1;
                    ej = edges_len(j);
                    ek = edges_len(mod(j,4)+1);
                    if (ei == 0 ) && (ej == 0)
                        %denerate face
                        p_i = ReservoirIntersection.rayPlaneIntersection(o,d,planes(i,:));
                        p_j = ReservoirIntersection.rayPlaneIntersection(o,d,planes(mod(j,4)+1,:));
                        new_coords(i,:) = p_i;
                        new_coords(i + 4,:) = p_i;
                        new_coords(j,:) = p_j;
                        new_coords(j + 4,:) = p_j;
                        break;
                    end
                end
                return
%                 view(162,14)
%                 ReservoirIntersection.showDegenerateCases(pts_old,planes,o,d);
%                 for i = 1:4
%                     ei = edges_len(i);
%                     j = mod(i,4)+1;
%                     ej = edges_len(j);
%                     if (ei == 0 ) 
%                         if (ej == 0)
%                             Drawer.drawPlan(planes(i,1:3)',planes(i,4),Drawer.m_blue);
%                             Drawer.drawPlan(planes(mod(j,4)+1,1:3)',planes(mod(j,4)+1,4), Drawer.m_blue);
%                         end
%                     else
%                         p_i = pts_old_order(i,:);
%                         p_j = pts_old_order(i+4,:);
%                         n = -cross(planes(i,1:3),planes(mod(i,4)+1,1:3));
%                         n = n / norm(n);
%                         fac = 1;
%                         p_i = p_i - fac*n;
%                         p_j = p_j + fac*n;
%                         line([p_i(1),p_j(1)],[p_i(2),p_j(2)],[p_i(3),p_j(3)],'linewidth',2.0,'color','magenta');
%                     end
%                 end
%                  for i = 1:4
%                      p_i = new_coords(i,:);
%                      p_j = new_coords(i+4,:);                     
%                      fac_i = 0.1;
%                      fac_j = 0.2;
%                      edge = -cross(planes(i,1:3),planes(mod(i,4)+1,1:3));;
%                      edge = edge/norm(edge);
%                      p_i_text = p_i - fac_i*edge;
%                      p_j_text = p_j + fac_j*edge;
%                      plot3(p_i(1),p_i(2),p_i(3),'o','MarkerSize',4,'MarkerFaceColor','green','Color','green');
%                      plot3(p_j(1),p_j(2),p_j(3),'o','MarkerSize',4,'MarkerFaceColor','green','Color','green');
%                      text(p_i_text(1),p_i_text(2),p_i_text(3),strcat('\boldmath$\hat{p}_{',num2str((i)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9,'Interpreter','latex');
%                      text(p_j_text(1),p_j_text(2),p_j_text(3),strcat('\boldmath$\hat{p}_{',num2str((i+4)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9,'Interpreter','latex');
%                  end
%                 exportgraphics(gca,strcat('deg_2','.png'),'Resolution',500);
%                 exportgraphics(gca,strcat('deg_2','.pdf'),'Resolution',500);
            end
%             if inverted ~= -1 && count == 1 
%                  view(27,69);
%                  ReservoirIntersection.showDegenerateCases(pts_old,planes,o,d);
%                  for j = 1:4
%                      if edges_len(j) == 0
%                          break
%                      end
%                  end
%                  Drawer.drawPlan(planes(j,1:3)',planes(j,4),Drawer.m_blue);
%                  Drawer.drawPlan(planes(mod(j,4)+1,1:3)',planes(mod(j,4)+1,4),Drawer.m_blue);
%                  %draw edges
%                  for i = 1:4
%                      p_i = pts_old_order(i,:);
%                      p_j = pts_old_order(i+4,:);
%                      n = -cross(planes(i,1:3),planes(mod(i,4)+1,1:3));
%                      n = n / norm(n);
%                      fac = 4;
%                      p_i = p_i - fac*n;
%                      p_j = p_j + fac*n;
%                      line([p_i(1),p_j(1)],[p_i(2),p_j(2)],[p_i(3),p_j(3)],'linewidth',1.5,'color','magenta');
%                  end
%                 for i = 1:4
%                     p_i = new_coords(i,:);
%                     p_j = new_coords(i+4,:);
%                     edge = p_j - p_i;
%                     edge = edge/norm(edge);
%                     fac_i = 0.2;
%                     fac_j = 0.5;
%                     p_i_text = p_i - fac_i*edge;
%                     p_j_text = p_j + fac_j*edge;
%                     plot3(p_i(1),p_i(2),p_i(3),'o','MarkerSize',4,'MarkerFaceColor','green','Color','green');
%                     plot3(p_j(1),p_j(2),p_j(3),'o','MarkerSize',4,'MarkerFaceColor','green','Color','green');
%                     text(p_i_text(1),p_i_text(2),p_i_text(3),strcat('\boldmath$\hat{p}_{',num2str((i)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9,'Interpreter','latex');
%                     text(p_j_text(1),p_j_text(2),p_j_text(3),strcat('\boldmath$\hat{p}_{',num2str((i+4)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9,'Interpreter','latex');
%                 end
%                  l = light('Position',[1 1 5.0],'Style','infinite');
%                  exportgraphics(gca,strcat('deg_1','.png'),'Resolution',500);
%                  exportgraphics(gca,strcat('deg_1','.pdf'),'Resolution',500);
%              end
%              if count == 3
%                  view(73,46);
%                  ReservoirIntersection.showDegenerateCases(pts_old,planes,o,d);
%                  Drawer.drawPlan(planes(1,1:3)',planes(1,4),Drawer.m_blue);
%                  Drawer.drawPlan(planes(2,1:3)',planes(2,4),Drawer.m_blue);
%                  Drawer.drawPlan(planes(3,1:3)',planes(3,4),Drawer.m_blue);
%                  Drawer.drawPlan(planes(4,1:3)',planes(4,4),Drawer.m_blue);
%                  %draw edges
%                  for i = 1:4
%                      p_i = pts_old_order(i,:);
%                      p_j = pts_old_order(i+4,:);
%                      n = cross(planes(i,1:3),planes(mod(i,4)+1,1:3));
%                      n = n / norm(n);
%                      fac = 3;
%                      p_i = p_i - fac*n;
%                      p_j = p_j + fac*n;
%                      line([p_i(1),p_j(1)],[p_i(2),p_j(2)],[p_i(3),p_j(3)],'linewidth',2.0,'color','magenta');
%                  end
%                  for i = 1:4
%                      p_i = new_coords(i,:);
%                      p_j = new_coords(i+4,:);
%                      edge = p_j - p_i;
%                      edge = edge/norm(edge);
%                      fac_i = 0.2;
%                      fac_j = 0.4;
%                      p_i_text = p_i - fac_i*edge;
%                      p_j_text = p_j + fac_j*edge;
%                      plot3(p_i(1),p_i(2),p_i(3),'o','MarkerSize',4,'MarkerFaceColor','green','Color','green');
%                      plot3(p_j(1),p_j(2),p_j(3),'o','MarkerSize',4,'MarkerFaceColor','green','Color','green');
%                      text(p_i_text(1),p_i_text(2),p_i_text(3),strcat('\boldmath$\hat{p}_{',num2str((i)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10,'Interpreter','latex');
%                      text(p_j_text(1),p_j_text(2),p_j_text(3),strcat('\boldmath$\hat{p}_{',num2str((i+4)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10,'Interpreter','latex');
%                  end
%                  %l = light('Position',[1 1 5.0],'Style','infinite');
%                  light('Position',[1 1 5.0],'Style','infinite');
%                  exportgraphics(gca,strcat('deg_3','.png'),'Resolution',500);
%                  exportgraphics(gca,strcat('deg_3','.pdf'),'Resolution',500);
%              end
%             if  (count == 0)
%                 view(14,32);
%                 ReservoirIntersection.showDegenerateCases(pts_old,planes,o,d);
%                 %draw edges
%                 for i = 1:4
%                     p_i = pts_old(i,:);
%                     p_j = pts_old(i+4,:);
%                     n = p_j - p_i;
%                     n = n / norm(n);
%                     fac = 4;
%                     p_i = p_i - fac*n;
%                     p_j = p_j + fac*n;
%                     line([p_i(1),p_j(1)],[p_i(2),p_j(2)],[p_i(3),p_j(3)],'linewidth',3,'color','magenta');
%                 end
%                 for i = 1:4
%                     p_i = new_coords(i,:);
%                     p_j = new_coords(i+4,:);
%                     edge = p_j - p_i;
%                     edge = edge/norm(edge);
%                     fac_i = 0.4;
%                     fac_j = 1.2;
%                     p_i_text = p_i - fac_i*edge;
%                     p_j_text = p_j + fac_j*edge;
%                     plot3(p_i(1),p_i(2),p_i(3),'o','MarkerSize',8,'MarkerFaceColor','green','Color','green');
%                     plot3(p_j(1),p_j(2),p_j(3),'o','MarkerSize',8,'MarkerFaceColor','green','Color','green');
%                     text(p_i_text(1),p_i_text(2),p_i_text(3),strcat('\boldmath$\hat{p}_{',num2str((i)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Interpreter','latex');
%                     text(p_j_text(1),p_j_text(2),p_j_text(3),strcat('\boldmath$\hat{p}_{',num2str((i+4)),'}$'), 'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Interpreter','latex');
%                 end
%                 l = light('Position',[1 1 5.0],'Style','infinite');
%                 exportgraphics(gca,strcat('deg_0','.png'),'Resolution',500);
%                 exportgraphics(gca,strcat('deg_0','.pdf'),'Resolution',500); 
%             end
        end
        function children_cells = childrenCells(pts)
            top = [pts(1,:);pts(2,:);pts(4,:); pts(3,:)];
            bottom = [pts(5,:);pts(6,:);pts(8,:); pts(7,:)];
            normals = ReservoirIntersection.createNormals(top);
            new_pts_top = ReservoirIntersection.subdivide(top);
            new_pts_bottom = ReservoirIntersection.subdivide(bottom);
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
        function [ratio,volume_total] = errConcave(pts)
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
            for i = 1:6
                coords = planes{i};
                A = coords{1}';
                B = coords{2}';
                C = coords{3}';
                D = coords{4}';
                volume = (1/6)*dot(A-R,cross((B-D),(C-A))) + (1/12)*dot(C-A,cross((D-A),(B-A)));
                volume_total = volume_total + volume;
                volume_concave = volume_concave + abs(((1/12)*dot(C-A,cross((D-A),(B-A)))));
                %volume4 = (1/6)*dot(B-R,cross((B-D),(C-A))) - (1/12)*dot(B-A,cross((C-A), (D-A))) + (1/6)*dot(B-A,cross((C-A),(B-D)));
                %volume7 = (1/6)*dot(B-R,cross((B-D),(C-A))) + (1/12)*dot(D-B,cross((A-B), (C-D))) ;
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
        function showDegenerateCases(pts_old,planes,o,d)
            %view(112,47);
            Drawer.showSurface(pts_old, 'cyan');
            Drawer.drawPlan(planes(5,1:3)',planes(5,4),Drawer.m_red);
            Drawer.drawPlan(planes(6,1:3)',planes(6,4),Drawer.m_red);
            Drawer.drawLine(o,d);
            Drawer.setBB(pts_old);
            %Drawer.Shadow(reshape(pts_old',24,1)',200);
            
        end
    end

end