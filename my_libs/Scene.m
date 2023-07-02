classdef Scene < handle
    properties 
        filename
        figure
        angleInit_az = 48; % 136 148
        angleInit_el = 15;
        count_images = 0;
        im = {};
        angle = 2;
        m_red = [0.9176 0.2627 0.2078]
        m_blue = [0.2588 0.5216 0.9569]
        m_type;
        drawRobot
        m_margin = 0.1;
        m_l = light('Position',[1 1 5.0],'Style','infinite'); %light
    end
    methods
        function this = Scene(filename, type)
            this.m_type = type;
            fig = figure;
            this.figure = fig;
            this.figure.Position = [50 100 1024  600];
            axis equal;
            axis off
            set(gca,'Color','None');

            %             lighting gouraud;
            %             a1 = camlight('right');
            %             a2 = camlight('left');
            %             a3 = camlight('headlight');
            this.filename = filename;
            set(gca, 'XTick', [], 'YTick', []);
            axis equal
            axis off
%             axis vis3d;
%             camproj('persp');
%             T = get(gca,'tightinset');
%             set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
            %set(gcf,'color','white');
            %lighting gouraud;
            view(this.angleInit_az, this.angleInit_el);
            hold on           
            switch type                
                case "SO3"
                    view(0,90);
                    this.drawSphere(pi,type);
                    this.drawRobot = @this.point;
                case "SE2"        
                    view(0,90);
                    this.drawBoard(type);
                    this.drawRobot = @this.drawRobot2D;
                case "SE3"
                    view(30,30);
                    %this.drawBoard(type);
                    this.drawRobot = @this.drawRobot3D;
            end
        end

        function update(this)
            frame = getframe(this.figure);
            this.im{end+1} = frame2im(frame);
            this.angleInit_az = this.angleInit_az + 1;
            view(this.angleInit_az, this.angleInit_el);
        end
        function get(this)
            frame = getframe(this.figure);
            this.im{end+1} = frame2im(frame);
        end
        function save(this)
            for idx = 1:length(this.im)
                [A,map] = rgb2ind(this.im{idx},256);
                if idx == 1
                    imwrite(A,map,this.filename,'gif','LoopCount',Inf,'DelayTime',1.5);
                elseif (idx > 1) && (idx < length(this.im))
                    imwrite(A,map,this.filename,'gif','WriteMode','append','DelayTime',0.08);
                else
                    imwrite(A,map,this.filename,'gif','WriteMode','append','DelayTime',180);
                end
            end
        end
        function reset(this)
            this.im = {};
        end
        function exportFrame(~,filename)
            exportgraphics(gca,strcat("presentation/",filename,'.jpeg'),'Resolution',300);
        end
        function plotP(this, P)
            plot3(P(:,1),P(:,2),P(:,3),'o','markerfacecolor',this.m_red,'markersize',7);
        end
        function plotVec(this, P,V,color)
            if nargin < 4
                color = 'black';
            end
            quiver3(P(:,1),P(:,2),P(:,3),V(:,1),V(:,2),V(:,3),0.5,'color',color,'linewidth',1);
        end
        function plotLine(this,Pi,Pj,color)
            if nargin < 4
                color = 'black';
            end
            line([Pi(1),Pj(1)],[Pi(2),Pj(2)],[Pi(3),Pj(3)],'color',color,'linewidth',1);
        end
        function h = plotGeodesic(this,pts,color,style)
            if nargin < 3
                color = 'black';
                style = '-';
            end
            if nargin < 4
                style = '-';
            end
            pt_0 = pts(1,:);
            pt_end = pts(end,:);            
            h = plot3(pts(:,1),pts(:,2),pts(:,3),'linewidth',2,'color',color,'linestyle',style);
            h(end+1) = plot3(pt_0(:,1),pt_0(:,2),pt_0(:,3),'o','markerfacecolor','#58C4DD','markersize',7);
            h(end+1) = plot3(pt_end(:,1),pt_end(:,2),pt_end(:,3),'o','markerfacecolor','#FC6255','markersize',7);
        end
        function h = drawPath(this,poses)             
             dim = poses(1).getDim();
             n = length(poses);
             h = [];
             for i = 1:n-1
                 pose = poses(i);
                 h1 = this.drawRobot(pose.m_data);
                 v = pose.logMap(pose.m_data, poses(i+1).m_data);
                 pts = pose.geodesic(v,pose.m_data);
                 t = pts(:,1:dim,end);
                 if (dim == 2)
                     t = [t,0*t(:,1)];
                 end
                 h2 = plot3(t(:,1),t(:,2),t(:,3),'.','color','black','markersize',10);
                 h = [h;h1;h2];
             end
             h1 = this.drawRobot(poses(i+1).m_data);
             h = [h;h1];
         end
         function bb = setBB(this,coords)
            margin = this.m_margin;
            x_min = min(coords(1,:)) - margin;
            x_max = max(coords(1,:)) + margin;
            y_min = min(coords(2,:)) - margin;
            y_max = max(coords(2,:)) + margin;
            z_min = min(coords(3,:)) - margin;
            z_max = max(coords(3,:)) + margin;
            bb = [x_min, x_max, y_min, y_max, z_min, z_max];
            axis(bb);
         end
         function h = drawRobot3D(this,T,color)
             if nargin == 2
                 color = "#699C52";
             end
             h = this.drawCube(T,color,0.7,this.m_l);
         end
    end
    methods (Static)   
        function removeMargins(ax)
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end
        function h = point(P)
            plot3(P(:,1),P(:,2),P(:,3),'o','markerfacecolor','#FC6255','markersize',7);
            %Scene.label(P);
        end
        function h = label(P,symbol)
            if nargin < 2
                symbol = "y";
            end
            h = [];
            for i = 1:size(P,1)
                name = strcat("$",symbol,"_",num2str(i-1),"$");
                h(end+1) = text(P(i,1),P(i,2),P(i,3),name,'horizontalalignment','right','verticalalignment','top','FontSize',15,'interpreter','latex');
            end
        end
        function h = descri(P,name,pose)
            if nargin < 3
                pose = 'top';
            end
            h = text(P(1,1),P(1,2),P(1,3),name,'horizontalalignment','left','verticalalignment',pose,'FontSize',10,'interpreter','latex');
        end
        function drawSphere(radius,name)
            if nargin == 0                
                radius = pi;
                name = 'M';
            end
            if nargin == 1
                name = 'M';
            end
            [X,Y,Z] = sphere(150);
            s = surf(radius*X,radius*Y,radius*Z,'EdgeAlpha',0.0,'FaceAlpha',0.3);
            colormap('gray');
            u = [1,1,1];
            %plotp(u/norm(u),0)
            %plotp(-u/norm(u),0)
            %drawLine(-u,u);
            txt_pose = 2.0*radius*u/norm(u);
            text(txt_pose(1),txt_pose(2),txt_pose(3),name,'horizontalalignment','right','verticalalignment','top','FontSize',15,'interpreter','latex');
            [X,Y,Z] = sphere(30);
            x = radius*X(:);
            y = radius*Y(:);
            z = radius*Z(:);
            conec = convhull(x,y,z);
            Scene.Shadow(conec,x,y,z);          
            %exportgraphics(gca,strcat('sphere','.jpeg'),'Resolution',333);
        end
%         function h = drawRobot2D(poses)
%             r = 0.25;
%             x = poses(:,1:3);
%             d = r*[cos(poses(:,4)),sin(poses(:,4)),0*poses(:,4)];
%             h = plot3(x(:,1),x(:,2),x(:,3),'o','MarkerFaceColor',[0.2588 0.5216 0.9569],'markersize',5,'markeredgecolor','black');
%             h(end+1) = quiver3(x(:,1),x(:,2),x(:,3),d(:,1),d(:,2),d(:,3),'linewidth',1,'color','black','autoscale','off','ShowArrowHead','off');
%             n = size(poses,1);
%             theta = linspace(0,2*pi,20)';
%         end
         function h = drawRobot2D(pose,color)
            if nargin == 1
                color = "#699C52";
            end
            %x = poses(:,1:3);
            %d = r*[cos(poses(:,4)),sin(poses(:,4)),0*poses(:,4)];
            %h = plot3(x(:,1),x(:,2),x(:,3),'o','MarkerFaceColor',[0.2588 0.5216 0.9569],'markersize',5,'markeredgecolor','black');
            %h(end+1) = quiver3(x(:,1),x(:,2),x(:,3),d(:,1),d(:,2),d(:,3),'linewidth',1,'color','black','autoscale','off','ShowArrowHead','off');
            %n = size(poses,1);            
            t = pose(1:2,3);
            r = 0.25;
            rot_matrix = pose(1:2,1:2);            
            thetas = linspace(0,2*pi,9);
            r = r*[cos(thetas);sin(thetas)];            
            r = rot_matrix*r + t;
            h = fill(r(1,:),r(2,:),color,'FaceAlpha',0.36);
            h(end+1) = line([t(1),r(1,3)],[t(2),r(2,3)],'color','black','linewidth',1); 
         end

         function h = drawLandmarks(landmarks,color)
             [m,n] = size(landmarks);
             if (m == 2)
                 landmarks = [landmarks; zeros(1,n);];
             end
             if nargin == 1
                 color = [0.2588 0.5216 0.9569];
             end
             h = plot3(landmarks(1,:),landmarks(2,:),landmarks(3,:),'diamond','MarkerSize',10,'markerfacecolor',color);
%              r = 0.25;
%              thetas = linspace(0,2*pi,5);
%              square = r*[cos(thetas);sin(thetas)];
%              h = [];
%              for i = 1:size(landmarks,2)
%                  r_i = square + landmarks(:,i);
%                  h(end+1) = fill(r_i(1,:),r_i(2,:),color,'FaceAlpha',0.7);
%              end
         end

         function drawBoard(name,min_z)            
            if nargin == 1
                min_z = -0.1;
            end
            % checkboard texture
            ch = repmat(1-0.2*xor((mod(repmat(0:128-1,128,1),8*2)>7), ...
            (mod(repmat((0:128-1)',1,128),8*2)>7)),[1 1 3])*0.5 + 0.5;
            x = linspace(-15,15,2);
            [X,Y] = meshgrid(x,x);
            Z = 0*X + min_z;
            sc = surf(X,Y,Z, ...
             'CData',ch,'FaceColor','texturemap', ...
             'SpecularStrength',0, 'DiffuseStrength',0, 'AmbientStrength',1);
            x = X(:);
            y = Y(:);
            z = Z(:);
            x(end+1) = 0;
            y(end+1) = 0;
            z(end+1) = 0;
            txt_pose = 2.5*[0,1,1] + 2*[0,0,max(z)];       
            %ext(txt_pose(1),txt_pose(2),txt_pose(3),name,'horizontalalignment','right','verticalalignment','top','FontSize',15,'interpreter','latex');
            %conec = delaunayTriangulation(x,y);
            %Scene.Shadow(conec.ConnectivityList,x,y,z);             
        end
        function tsh = Shadow(conec,x,y,z,l)
            nudge = 0.5;            
            %l = light('Position',[1 4 5.0],'Style','infinite');
            background_color = get(gcf,'Color');
            z_min = min(min(z));
            ground = [0,0,-1,z_min - nudge];
            %ground = [0,0,-1,1e3*-3.2680 + 1];
            light_pos = [l.Position strcmp(l.Style,'local')];
            d = ground * light_pos';
            shadow_mat = d*eye(4) - light_pos'*ground;
            U = [];
            h = height(x);
            U = [x,y,ones(h,1)*ground(4)];
            % for i = 1:size(p,1)
            %     V =  reshape(cells(i,:),3,8)';
            %     U_i = [V ones(size(V,1),1)]*shadow_mat';
            %     U_i = bsxfun(@rdivide,U_i(:,1:3),U_i(:,4));
            %     U(end+1:end+3,1:3) = U_i;
            % end
            color = 0.8*[1,1,1];
            tsh = trisurf(conec,U(:,1),U(:,2),U(:,3), ...
                'FaceColor',color, ...
                'DiffuseStrength',0,'SpecularStrength',0, ...
                'AmbientStrength',1, ...
                'EdgeColor','none');
            D = sum(bsxfun(@times,U(:,1:2),l.Position(1:2)),2);
            if numel(D)>1
                D = Scene.matrixnormalize(D);
            end
            %D = 1.0-D;
            C = bsxfun(@plus,color,bsxfun(@times,D,background_color-color));
            tsh.FaceVertexCData = C;
            tsh.FaceColor = 'interp';
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
        function h = drawCube(T,color,alpha,l)
            if nargin == 1
                color = "#699C52";
                alpha = 0.7;
            end
            if nargin == 2
                alpha = 0.7;
            end
            fac = 0.01/sqrt(3);
            p_max = fac * [1,1,1];            
            p_min = -p_max;          

            p5  = p_min;
            p6 = [p_max(1),p_min(2),p_min(3)];
            p8 = [p_max(1),p_max(2),p_min(3)];
            p7 = [p_min(1),p_max(2),p_min(3)];
            p1 = [p5(1),p5(2),p_max(3)];
            p2 = [p6(1),p6(2),p_max(3)];
            p4 = [p8(1),p8(2),p_max(3)];
            p3 = [p7(1),p7(2),p_max(3)];
            pts = [p1;p2;p3;p4;p5;p6;p7;p8];

            R = T(1:3,1:3);
            t = T(1:3,4)';
            pts = pts * R + t;
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
            h = line([pts(edges(:,1),1)';pts(edges(:,2),1)'],[pts(edges(:,1),2)';pts(edges(:,2),2)'],[pts(edges(:,1),3)';pts(edges(:,2),3)'],'Color',color,'linewidth',2.5);
            h = [h ; trisurf(faces,pts(:,1),pts(:,2),pts(:,3),'EdgeAlpha',0,'FaceColor',color,'FaceAlpha',alpha)];
            x = pts(:,1);
            y = pts(:,2);
            conec = delaunayTriangulation(x,y);
            Scene.Shadow(conec.ConnectivityList,x,y,0*x,l);             
        end
        function h = drawAxis(T)
            fac = 0.05;
            R = fac*T(1:3,1:3);
            o = T(1:3,4);
            h = quiver3(o(1),o(2),o(3),R(1,1),R(1,2),R(1,3),'linewidth',1,'color','red');
            h = [h, quiver3(o(1),o(2),o(3),R(2,1),R(2,2),R(2,3),'linewidth',1,'color','green')];
            h = [h, quiver3(o(1),o(2),o(3),R(3,1),R(3,2),R(3,3),'linewidth',1,'color','blue')];
        end
    end
end