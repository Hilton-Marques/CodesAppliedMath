classdef Drawer < handle
    properties (Constant)
        m_margin = 2;
        m_red = [0.9176 0.2627 0.2078];
        m_blue = [0.2588 0.5216 0.9569];
    end
    properties
        m_figure
        m_cell_intersection_obj CellIntersection
        m_reservoir_intersection_obj ReservoirIntersection
    end
    methods
        function this = Drawer(reservoir_intersection_obj)
            if nargin > 0
                %this.m_cell_intersection_obj = cell_intersection_obj;
                this.m_cell_intersection_obj = reservoir_intersection_obj;
                %this.m_reservoir_intersection_obj = reservoir_intersection_obj;
            end
            this.reset();
            %this.demoCells();
        end
        function demoCells(this)
            pts = [[-1,-1,1.0];
                [1,-1,2.0];
                [-1,1,2.0];
                [1,1,1];
                [-1,-1,-2.5];
                [1,-1,-1];
                [-1,1,-1];
                [1,1,-2.5]];
            this.showSurface(pts, this.m_red);
            axis vis3d;
            %camproj('persp');
            T = get(gca,'tightinset');
            pt_shadow = reshape(pts',24,1)';
            this.Shadow(pt_shadow,-1.5); %to pyramid
            %this.Shadow(pt_shadow,1.0);
            %view(162,14);
            view(-76,11); %to pyramid
            this.TitleExport(".pdf");
        end
        function showGrid(this,ids,alpha)
            n = size(ids,1);
            U = zeros(n,24);
            for i = 1:n
                bb = this.m_reservoir_intersection_obj.m_bb.GetBinBB(ids(i,:));
                pts = this.showBB(bb,alpha);
                U(i,:) = reshape(pts',24,1)';
            end
            this.Shadow(U,20);
        end
        function showModel(this, model, color)
            n = 30000;
            this.showCellss(model(1:n,:),color,0.5);
            axis vis3d;
            %camproj('persp');
            T = get(gca,'tightinset');
            set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
            n = 1000;
            this.Shadow(model(1:n,:));
            %this.TitleExport(".pdf");
        end

        function showParentChildren(this, parent, children)
            this.showSurface(parent, this.m_red);
            for i = 1:size(children,2)
                this.showSurface(children{i}, this.m_blue);
            end
        end
        
        function showCellss(this,cells,color,alpha)
            for i = 1:size(cells,1)
                coords = reshape(cells(i,:),3,8)';
                this.drawTopo(coords,color,alpha);
                %this.showSurface(coords, this.m_red);
            end
        end
        function showCells(this)
            title('cell A and cell B','Interpreter','latex')
            this.showSurface(this.m_cell_intersection_obj.m_cell_A, this.m_red);
            this.showSurface(this.m_cell_intersection_obj.m_cell_B, this.m_blue);
            %this.export();
            this.reset();
        end
        function showConvexApproximation(this)
            title('convex approximation','Interpreter','latex')
            %this.showSurface(this.m_cell_intersection_obj.m_cell_A, this.m_red,false);
            this.showSurface(this.m_cell_intersection_obj.m_cell_B, this.m_blue,false);

            %this.plotPlanes(this.m_cell_intersection_obj.m_planes_A(1:2,:),this.m_red);
            this.plotPlanes(this.m_cell_intersection_obj.m_planes_B(1:2,:),this.m_blue);
%             this.plotPlanes(this.m_cell_intersection_obj.m_planes_B(3:4,:),this.m_red);
%             this.plotPlanes(this.m_cell_intersection_obj.m_planes_B(5:6,:),this.m_blue);

            for i = 1:size(this.m_cell_intersection_obj.m_convex_A.pts,1)
                pt_i = this.m_cell_intersection_obj.m_convex_A.pts(i,:);
                pt_j = this.m_cell_intersection_obj.m_convex_B.pts(i,:);
                plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',7,'MarkerFaceColor','black');
                plot3(pt_j(1),pt_j(2),pt_j(3),'o','MarkerSize',7,'MarkerFaceColor','black');
                %text(pt_i(1),pt_i(2),pt_i(3),num2str(i-1), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',10);
            end
            %this.export();
            this.reset();
        end
        function showGJK(this)
            %title('GJK solution','Interpreter','latex')
            this.m_cell_intersection_obj.m_convex_A.show(this.m_red);
            this.m_cell_intersection_obj.m_convex_B.show(this.m_blue);
            plot3(this.m_cell_intersection_obj.m_p(1),this.m_cell_intersection_obj.m_p(2), this.m_cell_intersection_obj.m_p(3), 'o','MarkerSize',11,'MarkerFaceColor','black');
            U(1,:) = reshape(this.m_cell_intersection_obj.m_cell_A',24,1)';
            U(2,:) = reshape(this.m_cell_intersection_obj.m_cell_B',24,1)';
            this.Shadow(U,1);
            %this.export();
            this.reset();
        end
        function showVolumes(this)
            title('convex approximation','Interpreter','latex')
            this.showSurface(this.m_cell_intersection_obj.m_cell_A, this.m_red,false);
            this.show
            Surface(this.m_cell_intersection_obj.m_cell_B, this.m_blue,false);
        end
        function setVolumeTitle(this,volume)
            str = strcat("$volume =",num2str(volume,'%.2f'),"$");
            %title(str,'Interpreter','latex')
            this.export();
        end
        function TitleExport(this,str)
            if iscell(str)
                str = str{1};
            end
            %title(str,'Interpreter','latex', 'Units', 'normalized', 'Position', [0.5, 0.9, 0]);
            this.export(str);
        end
        function reset(this)
            this.m_figure = figure;
            hold on
            axis off
            set(gcf,'color','white');
            lighting gouraud;
            camlight('headlight');
            view(-9,48);
            if ~isempty(this.m_cell_intersection_obj)
                this.setBB(this.m_cell_intersection_obj.getCoords());
            end
            %camproj('persp');
        end
        function pts = showBB(this,bb,alpha)
            p_min = bb(1,:);
            p_max = bb(2,:);
            p5  = p_min;
            p6 = [p_max(1),p_min(2),p_min(3)];
            p8 = [p_max(1),p_max(2),p_min(3)];
            p7 = [p_min(1),p_max(2),p_min(3)];

            p1 = [p5(1),p5(2),p_max(3)];
            p2 = [p6(1),p6(2),p_max(3)];
            p4 = [p8(1),p8(2),p_max(3)];
            p3 = [p7(1),p7(2),p_max(3)];

            pts = [p1;p2;p3;p4;p5;p6;p7;p8];
            this.(pts,0.75*[1,1,1],alpha);
        end
    end
    methods (Static)
        function Shadow(cells,nudge)
            if nargin == 2
                nudge = 2e-3;
            end
            %l = light('Position',[1 4 5.0],'Style','infinite');
            l = light('Position',[1 1 5.0],'Style','infinite');
            background_color = get(gcf,'Color');
            Faces = [[1,2,4];
                [4,3,1];%top
                [5,6,8];
                [8,7,5];%bottom
                [6,8,4];
                [4,2,6]; %right
                [1,3,7];
                [7,5,1];%left
                [8,4,3];
                [3,7,8];%back
                [5,6,2];
                [2,1,5]];%front
            z =cells(:,3:3:24);
            z_min = min(min(z));
            ground = [0,0,-1,z_min - nudge];
            %ground = [0,0,-1,1e3*-3.2680 + 1];
            light_pos = [l.Position strcmp(l.Style,'local')];
            d = ground * light_pos';
            shadow_mat = d*eye(4) - light_pos'*ground;
            U = [];
            conec = [];
            for i = 1:size(cells,1)
                V =  reshape(cells(i,:),3,8)';
                U_i = [V ones(size(V,1),1)]*shadow_mat';
                U_i = bsxfun(@rdivide,U_i(:,1:3),U_i(:,4));
                U(end+1:end+8,1:3) = U_i;
                conec(end+1:end+12,1:3) =  Faces + (i-1)*8;
            end
            color = 0.8*[1,1,1];
            tsh = trisurf(conec,U(:,1),U(:,2),U(:,3), ...
                'FaceColor',color, ...
                'DiffuseStrength',0,'SpecularStrength',0, ...
                'AmbientStrength',1, ...
                'EdgeColor','none');
            D = sum(bsxfun(@times,U(:,1:2),l.Position(1:2)),2);
            if numel(D)>1
                D = Drawer.matrixnormalize(D);
            end
            %D = 1.0-D;
            C = bsxfun(@plus,color,bsxfun(@times,D,background_color-color));
            tsh.FaceVertexCData = C;
            tsh.FaceColor = 'interp';
        end
        function setBB(coords)
            margin = Drawer.m_margin;
            x_min = min(coords(:,1)) - margin;
            x_max = max(coords(:,1)) + margin;
            y_min = min(coords(:,2)) - margin;
            y_max = max(coords(:,2)) + margin;
            z_min = min(coords(:,3)) - margin;
            z_max = max(coords(:,3)) + margin;
            bb = [x_min,x_max,y_min,y_max,z_min,z_max];
            axis(bb);
        end
        function showSurface(pts,color,flag)
            if nargin == 2
                flag = true;
            end
            n = 10;
            ksi = linspace(0,1,n);
            eta = linspace(0,1,n);
            [ksi,eta] = meshgrid(ksi,eta);
            planes = cell(6,1);
            planes{1} = {pts(1,:),pts(2,:),pts(4,:), pts(3,:)};
            planes{2} = {pts(5,:),pts(7,:),pts(8,:), pts(6,:)};
            planes{3} = {pts(1,:),pts(3,:),pts(7,:), pts(5,:)};
            planes{4} = {pts(2,:),pts(6,:),pts(8,:), pts(4,:)};
            planes{5} = {pts(1,:),pts(5,:),pts(6,:), pts(2,:)};
            planes{6} = {pts(3,:),pts(4,:),pts(8,:), pts(7,:)};
            ids = [0,1,3,2];
            if flag
                for i = 1:size(pts,1)
                    %for i = 1:4
                    pt_i = pts(i,:);
                    %plot3(pt_i(1),pt_i(2),pt_i(3),'o','MarkerSize',4,'MarkerFaceColor','black','Color','black');
                    if i <= 4
                        %text(pt_i(1),pt_i(2),pt_i(3),strcat('\boldmath$p_',num2str(ids(i)),'$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',12,'Interpreter','latex');
                    end
                end
            end
            k = [1,2,3;1,3,4;1,4,5;1,5,2];
            R = sum(pts,1)/8;
            %plot3(R(1),R(2),R(3),'o','MarkerSize',4,'MarkerFaceColor','black','Color','black');
            %text(R(1),R(2),R(3),strcat('\boldmath$r$'), 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',18,'Interpreter','latex');
            for i = 1:6%1,6
                coords = planes{i};
                p0 = coords{1}';
                p1 = coords{2}';
                p2 = coords{3}';
                p3 = coords{4}';
                x = p0(1) + ksi.*(p1(1) - p0(1)) + eta.*(p3(1) - p0(1)) + ksi.*eta.*( (p2(1)-p0(1)) - (p1(1)-p0(1)) - (p3(1) - p0(1)) );
                y = p0(2) + ksi.*(p1(2) - p0(2)) + eta.*(p3(2) - p0(2)) + ksi.*eta.*((p2(2)-p0(2)) - (p1(2)-p0(2)) - (p3(2) - p0(2)));
                z = p0(3) + ksi.*(p1(3) - p0(3)) + eta.*(p3(3) - p0(3)) + ksi.*eta.*((p2(3)-p0(3)) - (p1(3)-p0(3)) - (p3(3) - p0(3)));
                surf(x,y,z,'FaceAlpha',0.7,'EdgeColor','none','FaceColor', color);
                pyramid = [R;p0';p1';p2';p3'];
                %kk = convhull(pyramid(:,1),pyramid(:,2),pyramid(:,3));

                %trisurf(k,pyramid(:,1),pyramid(:,2),pyramid(:,3),'FaceColor',0.7*ones(1,3),'FaceAlpha',0.8);%0.8
                %trisurf(kk,pyramid(:,1),pyramid(:,2),pyramid(:,3),'FaceColor','cyan');
                %tetrahedra
                tetra = [p0';p1';p3';p1';p2';p3'];
                %trisurf([[1,2,3];[4,5,6]], tetra(:,1),tetra(:,2),tetra(:,3),'FaceColor',0.7*ones(1,3),'FaceAlpha',0.6);
                tetra = [p0';p1';p2';p3'];
                %trisurf([[1,2,4];[2,3,4];[1,3,4];[1,2,3]],tetra(:,1),tetra(:,2),tetra(:,3),'FaceColor',0.7*ones(1,3),'FaceAlpha',0.0);
                line([p0(1),p1(1)],[p0(2),p1(2)],[p0(3),p1(3)],'color','black');
                line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'color','black');
                line([p2(1),p3(1)],[p2(2),p3(2)],[p2(3),p3(3)],'color','black');
                line([p3(1),p0(1)],[p3(2),p0(2)],[p3(3),p0(3)],'color','black');
            end
        end

        function showFace(plan)
            n = 20;
             ksi = linspace(0,1,n);
            eta = linspace(0,1,n);
            [ksi,eta] = meshgrid(ksi,eta);
            coords = plan;
            p0 = coords(1,:)';
            p1 = coords(2,:)';
            p2 = coords(3,:)';
            p3 = coords(4,:)';
            x = p0(1) + ksi.*(p1(1) - p0(1)) + eta.*(p3(1) - p0(1)) + ksi.*eta.*( (p2(1)-p0(1)) - (p1(1)-p0(1)) - (p3(1) - p0(1)) );
            y = p0(2) + ksi.*(p1(2) - p0(2)) + eta.*(p3(2) - p0(2)) + ksi.*eta.*((p2(2)-p0(2)) - (p1(2)-p0(2)) - (p3(2) - p0(2)));
            z = p0(3) + ksi.*(p1(3) - p0(3)) + eta.*(p3(3) - p0(3)) + ksi.*eta.*((p2(3)-p0(3)) - (p1(3)-p0(3)) - (p3(3) - p0(3)));
            surf(x,y,z,'FaceAlpha',0.7,'EdgeColor','none','FaceColor', 'green');
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
           % h1 = fill3([xp1(1),yp1(1),xp2(1),yp2(1)],...
             %   [xp1(2),yp1(2),xp2(2),yp2(2)],...
             %   [xp1(3),yp1(3),xp2(3),yp2(3)],color);
            ids = [1,2,3;1,3,4];
            pts = [xp1';yp1';xp2';yp2'];
            h1  =trisurf(ids,pts(:,1),pts(:,2),pts(:,3),'EdgeColor','none','FaceColor',color,'FaceAlpha',0.8);
            set(h1, 'facealpha',0.8);
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
            line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'linewidth',2,'color','black');
        end
        function export(type)
            if nargin < 1
                type = '.png';
            end
            str = get(gca,'title').String;
            %view(35,41);
            view(56,23);
            exportgraphics(gca,strcat(str,type),'Resolution',2000);%'ContentType','vector'
        end
        
        function drawTopo(pts,color,alpha)
            if (nargin == 2)
                alpha = 0.6;
            end
            edges = [[1,2];
                [2,4];
                [4,3];
                [3,1];
                [5,6];
                [6,8];
                [8,7];
                [7,5];
                [5,1];
                [6,2];
                [8,4];
                [7,3]];

            faces = [[1,2,4];
                [4,3,1];%top
                [5,6,8];
                [8,7,5];%bottom
                [6,8,4];
                [4,2,6]; %right
                [1,3,7];
                [7,5,1];%left
                [8,4,3];
                [3,7,8];%back
                [5,6,2];
                [2,1,5]];%front
            line([pts(edges(:,1),1)';pts(edges(:,2),1)'],[pts(edges(:,1),2)';pts(edges(:,2),2)'],[pts(edges(:,1),3)';pts(edges(:,2),3)'],'Color',color,'linewidth',0.5);
            %tf = trisurf(faces,pts(:,1),pts(:,2),pts(:,3),'EdgeAlpha',0,'FaceColor',color,'FaceAlpha',alpha);
        end

        function [h,L,M,ground] = add_shadow(T,L,varargin)
            % ADD_SHADOW  Add a shadow for plotted mesh (t) according to light (l)
            %
            % h = add_shadow()  % Apply to all
            % h = add_shadow(T,L)
            % h = add_shadow(T,L,'ParameterName',ParameterValue,...)
            % h = add_shadow([],[],'ParameterName',ParameterValue,...) % apply to all
            %
            % Inputs:
            %   T  #T list of trisurf handles {[] --> find all trisurf in `gca`}
            %   L  #L list of lights {[] --> find all light in `gca`}
            %   Optional:
            %     'Ground'  ground plane equation {[0 0 -1 min(Z)]}
            %     'Nudge'  nudge the ground plane down a bit
            %     'Color' followed by 3-vector color {get(gcf,'Color')*0.9}
            %     'BackgroundColor' followed by 3-vector color {get(gcf,'Color')*0.9}
            %     'Fade'  followed by:
            %        'none' constant shadow color
            %        'local' fade darker away from contact with ground (ape a spotlight)
            %        {'infinite'} fade lighter away from contact ground (ape infinite
            %          light)
            % Outputs:
            %   h  #T*#L list of output shadow trisurf handles
            %   L  #L list of lights
            %   M  4 by 4 by #T*#L shadow projection matrices
            %
            % Example:
            %   t = tsurf(F,V,'EdgeColor','none',fphong,'SpecularStrength',0.1);
            %   l = [
            %      light('Position',[-10 -10 13],'Style','local');
            %      light('Position',[10 -10  13],'Style','local')];
            %   camproj('persp');
            %   axis equal;
            %   h = add_shadow(t,l);
            %   apply_ambient_occlusion(t);
            %

            % default values
            ground = [];
            nudge = 0;
            color = get(gcf,'Color')*0.9;
            background_color = get(gcf,'Color');
            fade = 'infinite';
            % Map of parameter names to variable names
            params_to_variables = containers.Map( ...
                {'Ground','Nudge','BackgroundColor','Color','Fade'}, ...
                {'ground','nudge','background_color','color','fade'});
            v = 1;
            while v <= numel(varargin)
                param_name = varargin{v};
                if isKey(params_to_variables,param_name)
                    assert(v+1<=numel(varargin));
                    v = v+1;
                    % Trick: use feval on anonymous function to use assignin to this workspace
                    feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
                else
                    error('Unsupported parameter: %s',varargin{v});
                end
                v=v+1;
            end

            if ~exist('T','var') || isempty(T)
                c = get(gca,'Children');
                T = c(arrayfun(@(x) ...
                    isa(x,'matlab.graphics.primitive.Patch') || ...
                    isa(x,'matlab.graphics.chart.primitive.Scatter') ...
                    ,c));
            end
            if iscell(T)
                T = [T{:}]';
            end
            T = T(:);

            ca = caxis;
            if ~exist('L','var') || isempty(L)
                c = get(gca,'Children');
                L = c(arrayfun(@(x) isa(x,'matlab.graphics.primitive.Light'),c));
                if isempty(L)
                    L = camlight;
                end
            end

            if isempty(ground)
                minZ = inf;
                for t = T'
                    switch class(t)
                        case 'matlab.graphics.primitive.Patch'
                            V = t.Vertices;
                        case 'matlab.graphics.chart.primitive.Scatter'
                            V = [t.XData;t.YData;t.ZData]';
                        otherwise
                            class(t)
                    end
                    minZ = min([V(:,3);minZ]);
                end
                ground = [0 0 -1 minZ];
            end
            ground(4) = ground(4)-nudge;

            h = {};
            % need to specify that there are 0 "tubes"
            M = zeros([0,0,0]);
            for t = T'
                switch class(t)
                    case 'matlab.graphics.primitive.Patch'
                        V = t.Vertices;
                    case 'matlab.graphics.chart.primitive.Scatter'
                        V = [t.XData;t.YData;t.ZData]';
                    otherwise
                        class(t)
                end
                for l = L'
                    % plane equation
                    % 0 = ax + by + cz + d
                    light_pos = [l.Position strcmp(l.Style,'local')];
                    d = ground * light_pos';
                    shadow_mat = d*eye(4) - light_pos'*ground;
                    U = [V ones(size(V,1),1)]*shadow_mat';
                    U = bsxfun(@rdivide,U(:,1:3),U(:,4));

                    hold on;
                    switch class(t)
                        case 'matlab.graphics.primitive.Patch'
                            tsh = trisurf(t.Faces,U(:,1),U(:,2),U(:,3), ...
                                'FaceColor',color, ...
                                'DiffuseStrength',0,'SpecularStrength',0, ...
                                'AmbientStrength',1, ...
                                'EdgeColor','none');
                        case 'matlab.graphics.chart.primitive.Scatter'
                            tsh = copyobj(t,t.Parent);
                            tsh.XData = U(:,1);
                            tsh.YData = U(:,2);
                            tsh.ZData = U(:,3);
                            tsh.MarkerFaceColor = color;
                            tsh.MarkerEdgeColor = color;
                    end
                    hold off;
                    %continue

                    switch fade
                        case {'local','infinite'}
                            D = sum(bsxfun(@times,U(:,1:2),l.Position(1:2)),2);
                            if numel(D)>1
                                D = Drawer.matrixnormalize(D);
                            end
                            switch fade
                                case 'infinite'
                                    D = 1.0-D;
                            end
                            C = bsxfun(@plus,color,bsxfun(@times,D,background_color-color));
                            switch class(t)
                                case 'matlab.graphics.primitive.Patch'
                                    tsh.FaceVertexCData = C;
                                    tsh.FaceColor = 'interp';
                                case 'matlab.graphics.chart.primitive.Scatter'
                                    tsh.MarkerEdgeColor = 'flat';
                                    tsh.MarkerFaceColor = 'flat';
                                    tsh.CData = C;
                            end
                            caxis(ca);
                    end
                    % wireframe
                    switch class(t)
                        case 'matlab.graphics.primitive.Patch'
                            if t.FaceAlpha == 0 || (ischar(t.FaceColor) & strcmp(t.FaceColor,'none'))
                                tsh.FaceAlpha = 0;
                                tsh.EdgeColor = color;
                            end
                    end

                    h = {h{:} tsh};
                    M(:,:,end+1) = shadow_mat;
                end
            end
        end
        function N = matrixnormalize(M)
            % MATRIXNORMALIZE Normalize matrix values to be between the range 0 and 1.
            % Current just works with matrices of type double.
            %
            % N = matrixnormalize(M)
            %
            % Inputs:
            %   M  original input matrix
            % Output:
            %   N  normalized matrix
            %
            % Example:
            %   imshow(matrixnormalize(im))
            %   % Equivalent to
            %   imshow(im,[])
            %
            switch class(M)
                case {'double','single'}
                    N = (M-min(M(:)))./(max(M(:))-min(M(:)));
                case 'logical'
                    N = M;
                case 'uint8'
                    N = (M-min(M(:)))*(255/double((max(M(:))-min(M(:)))));
                otherwise
                    error('Class not supported');
            end
        end
    end
end