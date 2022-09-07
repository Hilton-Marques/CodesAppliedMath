classdef Drawer < handle
    properties
        m_margin = 10;
        m_figure
        m_red = [0.9176 0.2627 0.2078];
        m_blue = [0.2588 0.5216 0.9569];
        m_cell_intersection_obj CellIntersection

    end
    methods
        function this = Drawer(cell_intersection_obj)
            this.m_cell_intersection_obj = cell_intersection_obj;
            this.reset();
        end
        function setBB(this, coords)
            margin = this.m_margin;
            x_min = min(coords(:,1)) - margin;
            x_max = max(coords(:,1)) + margin;
            y_min = min(coords(:,2)) - margin;
            y_max = max(coords(:,2)) + margin;
            z_min = min(coords(:,3)) - margin;
            z_max = max(coords(:,3)) + margin;
            bb = [x_min,x_max,y_min,y_max,z_min,z_max];
            axis(bb);
        end
        function showParentChildren(this, parent, children)
            this.showSurface(parent, this.m_red);
            for i = 1:size(children,2)
                this.showSurface(children{i}, this.m_blue);
            end
        end
        function showCells(this)
            title('cell A and cell B','Interpreter','latex')
            this.showSurface(this.m_cell_intersection_obj.m_cell_A, this.m_red);
            this.showSurface(this.m_cell_intersection_obj.m_cell_B, this.m_blue);
            this.export();
            this.reset();
        end
        function showConvexApproximation(this)
            title('convex approximation','Interpreter','latex')
            this.showSurface(this.m_cell_intersection_obj.m_cell_A, this.m_red,false);
            this.showSurface(this.m_cell_intersection_obj.m_cell_B, this.m_blue,false);
            
            this.plotPlanes(this.m_cell_intersection_obj.m_planes_A(1:2,:),this.m_red);
            this.plotPlanes(this.m_cell_intersection_obj.m_planes_B(1:2,:),this.m_blue);

            for i = 1:size(this.m_cell_intersection_obj.m_convex_A,1)
                pt_i = this.m_cell_intersection_obj.m_convex_A(i,:);
                pt_j = this.m_cell_intersection_obj.m_convex_B(i,:);
                plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',7,'MarkerFaceColor','black');
                plot3(pt_j(1),pt_j(2),pt_j(3),'o','MarkerSize',7,'MarkerFaceColor','black');
                %text(pt_i(1),pt_i(2),pt_i(3),num2str(i-1), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
            end
            this.export();
            this.reset();
        end
        function showGJK(this)
            title('GJK solution','Interpreter','latex')
            this.m_cell_intersection_obj.m_convex_A.show(this.m_red);
            this.m_cell_intersection_obj.m_convex_B.show(this.m_blue);
            plot3(this.m_cell_intersection_obj.m_p(1),this.m_cell_intersection_obj.m_p(2), this.m_cell_intersection_obj.m_p(3), 'o','MarkerSize',10,'MarkerFaceColor','black');
            this.export();
            this.reset();
        end
        function showVolumes(this)
            title('convex approximation','Interpreter','latex')
            this.showSurface(this.m_cell_intersection_obj.m_cell_A, this.m_red,false);
            this.showSurface(this.m_cell_intersection_obj.m_cell_B, this.m_blue,false);
        end
        function setVolumeTitle(this,volume)
            str = strcat("$volume =",num2str(volume,'%.2f'),"$");
            title(str,'Interpreter','latex')
            this.export();
        end
        function reset(this)
            this.m_figure = figure;
            hold on
            axis off
            set(gcf,'color','white');
            lighting gouraud;
            camlight('headlight');
            view(-9,48);
            this.setBB(this.m_cell_intersection_obj.getCoords());
        end
    end
    methods (Static)
        function showSurface(pts,color,flag)
            if nargin == 2
                flag = true;
            end
            ksi = linspace(0,1,10);
            eta = linspace(0,1,10);
            [ksi,eta] = meshgrid(ksi,eta);
            planes = cell(6,1);
            planes{1} = {pts(1,:),pts(2,:),pts(4,:), pts(3,:)};
            planes{2} = {pts(5,:),pts(7,:),pts(8,:), pts(6,:)};
            planes{3} = {pts(1,:),pts(3,:),pts(7,:), pts(5,:)};
            planes{4} = {pts(2,:),pts(6,:),pts(8,:), pts(4,:)};
            planes{5} = {pts(1,:),pts(5,:),pts(6,:), pts(2,:)};
            planes{6} = {pts(3,:),pts(4,:),pts(8,:), pts(7,:)};
            if flag
                for i = 1:size(pts,1)
                    pt_i = pts(i,:);
                    plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',5,'MarkerFaceColor','r');
                    %text(pt_i(1),pt_i(2),pt_i(3),num2str(i-1), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
                end
            end
            k = [1,2,3;1,3,4;1,4,5;1,5,2];
            R = sum(pts,1)/8;
            plot3(R(1),R(2),R(3),'o','MarkerSize',4,'MarkerFaceColor','b');
            for i = 1:6
                coords = planes{i};
                p0 = coords{1}';
                p1 = coords{2}';
                p2 = coords{3}';
                p3 = coords{4}';

                x = p0(1) + ksi.*(p1(1) - p0(1)) + eta.*(p3(1) - p0(1)) + ksi.*eta.*( (p2(1)-p0(1)) - (p1(1)-p0(1)) - (p3(1) - p0(1)) );
                y = p0(2) + ksi.*(p1(2) - p0(2)) + eta.*(p3(2) - p0(2)) + ksi.*eta.*((p2(2)-p0(2)) - (p1(2)-p0(2)) - (p3(2) - p0(2)));
                z = p0(3) + ksi.*(p1(3) - p0(3)) + eta.*(p3(3) - p0(3)) + ksi.*eta.*((p2(3)-p0(3)) - (p1(3)-p0(3)) - (p3(3) - p0(3)));
                surf(x,y,z,'FaceAlpha',0.5,'EdgeAlpha',0.0,'FaceColor', color);
                pyramid = [R;p0';p1';p2';p3'];
                trisurf(k,pyramid(:,1),pyramid(:,2),pyramid(:,3),'FaceColor',color,'FaceAlpha',0.0);
            end
        end
        function plotPlanes(planes,color)
            for i = 1:size(planes,1)
                n = planes(i,1:3);
                d = planes(i,4);
                Drawer.drawPlan(n',d,color);
            end
        end
        function h1 = drawPlan(n,d,color,A,xf,xo)
            if nargin == 1
                A = 20000;
                xo = [0;0;0];
                color = [1,0,0];
                d = 0;
                xf = -n*d;
            elseif nargin == 2
                xo = [0;0;0];
                color = [1,0,0];
                A = 20000;
                xf = -n*d;
            elseif nargin == 3
                xo = [0;0;0];
                A = 2000000;
                xf = -n*d;
            elseif nargin == 5
                xo = [0;0;0];
            end
            [x,y] = Drawer.findTriedro(n);
            xp1 = xo + A*x + xf ;
            yp1 = xo + A*0.1*y + xf ;
            xp2 = xo - A*x + xf;
            yp2 = xo - A*0.1*y + xf;
            h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
                [xp1(2),yp1(2),xp2(2),yp2(2)],...
                [xp1(3),yp1(3),xp2(3),yp2(3)],color);
            set(h1, 'facealpha',0.5);
        end
        function [x,y] = findTriedro(z)
            z = z/norm(z);
            uTemp = [z(3);z(1);-z(2)];
            uTemp = uTemp/norm(uTemp);
            if (dot(z,uTemp) == 1)
                uTemp = uTemp([2,1,3]);
            end
            x = cross(uTemp,z);
            y = cross(z,x);
        end
        function plotPt(P,color)
            if nargin == 1
                color = 'black';
            end
            plot3(P(1),P(2),P(3),'o','MarkerFaceColor',color,'MarkerSize',5);
        end
        function drawLine(o,d)
            lam = 10000;
            p1 = o + lam*d;
            p2 = o - lam*d;
            line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'linewidth',3,'color','black');
        end
        function export()
            str = get(gca,'title').String;
            exportgraphics(gca,strcat(str,'.png'),'Resolution',333);
        end
    end
end