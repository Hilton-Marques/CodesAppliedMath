classdef Solid < handle
    properties
        points;
        edges;
        heds;
        elementsMother;
        elements = Face.empty;
        nL;
        nHeds;
        nPts;
        nEdges;
        nEl;
        elementsLevel;
        pointsBd;
        centroide;
    end
    methods
        function this = Solid(points,edges,heds,elements,nL)
            this.points = points;
            centroide = zeros(3,1);
            for pt = points
                centroide = centroide + pt.coord;
            end
            this.centroide = centroide/length(points);
            this.edges = edges;
            this.heds = heds;
            this.elementsMother = elements;
            this.nL = nL;
            this.nHeds = length(this.heds);
            this.nPts = length(this.points);
            this.nEdges = length(this.edges);
            this.nEl = length(elements);
            this.elementsLevel = cell(nL,1);
            this.start();
            
            this.initializeEl();
            this.initializePoints();
        end
        function start(this)
            this.elementsLevel{1} = this.elementsMother;
            for i = 2:this.nL
                Elements_m = this.elementsLevel{i-1};
                %update vertices
                for pt = this.points
                    star = pt.getStar(this);
                    n = length(star);
                    star_points = this.points(star);
                    %beta = 3/(8*n);
                    if n == 3
                        beta = (3/16);
                    else
                        %beta = 1/n*( 5/8 - (3/8 + 1/4*cos(2*pi/n))^2 );
                        beta = 3/(8*n);
                    end
                    pt.new_coord = (1 - n*beta)*pt.coord + beta*sum([star_points.coord],2);
                    %pt.coord = pt.new_coord;
                end
                npts = length(this.points);
                nElNv = length(Elements_m);
                newElements(4*nElNv) = Face;
                count = 1;
                for j = 1:nElNv
                    el = Elements_m(j);
                    localHeds = Hed.empty;
                    localPoints = Point.empty;
                    for k = 1:length(el.heds)
                        hed_m = el.heds(k);
                        edge_m = this.edges(hed_m.edgeId);
                        if ~edge_m.isSplited
                            p3 = this.points(this.heds(this.edges(hed_m.edgeId).getTwin(hed_m.id).heNext).inc(2)).coord;
                            p2 = this.points(this.heds(hed_m.heNext).inc(2)).coord;
                            p0 = this.points(hed_m.inc(1)).coord;
                            p1 = this.points(hed_m.inc(2)).coord;
                            coord = (p3 + p2 + 3*p0 + 3*p1)/8;
                            this.nPts = this.nPts + 1;
                            pt = Point(coord,this.nPts);
                            pts_heds = [this.points(hed_m.inc(1)),pt];
                            for l = 1:2 % number of subdivisions
                                this.nEdges = this.nEdges + 1;
                                this.nHeds = this.nHeds + 1;
                                inc1 = [hed_m.inc(l),this.nPts];
                                inc2 = [this.nPts,hed_m.inc(l)];
                                incs = [inc1;inc2];
                                newHed = Hed(incs(l,:),this.nHeds,this.nEdges);
                                this.nHeds = this.nHeds + 1;
                                newHedTwin = Hed(incs(mod(l,2)+1,:),this.nHeds,this.nEdges);
                                newEdge = Edge(newHed,newHedTwin,this.nEdges);
                                localHeds(end+1) = newHed;
                                this.heds(end+1) = newHed;
                                this.heds(end+1) = newHedTwin;
                                this.edges(end+1) = newEdge;
                                edge_m.childrenId(l) = this.nEdges;
                                pts_heds(l).hed0 = newHed;
                            end
                            %pt.hed0 = newHed;
                            edge_m.isSplited = true;
                            edge_m.mid = pt;
                            this.points(end+1) = pt;
                            localPoints(end+1) = pt;
                        else
                            localPoints(end+1) = edge_m.mid;
                            for l = 2:-1:1 % the this.heds reverse the order
                                edge_new = this.edges(edge_m.childrenId(l));
                                localHeds(end+1) = edge_new.hed2;
                            end
                        end
                    end
                    internalHeds(3) = Hed;
                    for k = 1:3
                        this.nEdges = this.nEdges + 1;
                        this.nHeds = this.nHeds + 1;
                        p2 = localPoints(mod(k,3)+1).id;
                        p1 = localPoints(k).id;
                        newHed = Hed([p1,p2],this.nHeds,this.nEdges);
                        this.nHeds = this.nHeds + 1;
                        newHedTwin = Hed([p2,p1],this.nHeds,this.nEdges);
                        internalHeds(k) = newHed;
                        this.heds(end+1) = newHed;
                        this.heds(end+1) = newHedTwin;
                        newEdge = Edge(newHed,newHedTwin,this.nEdges);
                        this.edges(end+1) = newEdge;
                    end
                    %reorder
                    internalHeds = internalHeds([3,1,2]) ;
                    % Update hNext for border triangles
                    c1 = 1;
                    c2 = 5;
                    for l = 1:2:5
                        this.nEl = this.nEl + 1;
                        firstId = localHeds(l).id;
                        idMidle = this.edges(internalHeds(c1).edgeId).hed2.id;
                        localHeds(l).heNext = idMidle;
                        localHeds(l).elId = this.nEl;
                        idBefore = mod(6 + l - 2,6)+1;
                        this.heds(idMidle).heNext = localHeds(idBefore).id;
                        this.heds(idMidle).elId = this.nEl;
                        localHeds(mod(c2,6)+1).heNext = firstId;
                        localHeds(mod(c2,6)+1).elId = this.nEl;
                        c1 = c1 + 1;
                        c2 = c2 + 2;
                    end
                    % Update hNext for internal triangle
                    this.nEl = this.nEl + 1;
                    for l = 1:3
                        internalHeds(l).heNext = internalHeds(mod(l,3)+1).id;
                        internalHeds(l).elId = this.nEl;
                    end
                    %Create new Elements
                    elementsChildren(4) = Face;
                    elementsChildren(1) = Face(localHeds(1).id,this.heds,el.q,el.typeEl);
                    elementsChildren(1).father = el;
                    elementsChildren(2) = Face(localHeds(3).id,this.heds,el.q,el.typeEl);
                    elementsChildren(2).father = el;
                    elementsChildren(3) = Face(localHeds(5).id,this.heds,el.q,el.typeEl);
                    elementsChildren(3).father = el;
                    elementsChildren(4) = Face(internalHeds(1).id,this.heds,el.q,el.typeEl);
                    elementsChildren(4).father = el;
                    newElements(count:count+3) = elementsChildren;
                    el.children = elementsChildren;
                    count = count + 4;
                end
                % Create new nodes
                for el = newElements
                    newPoints = el.refine(this.edges,this.points);
                    % As new this.points was created there is a need to updtate this.nPts
                    this.points(end+1:end+length(newPoints)) = newPoints;
                    this.nPts = length(this.points);
                end
                this.elementsLevel{i} = newElements;
                for v = 1:npts
                    this.points(v).coord = this.points(v).new_coord;
                end
                
            end
            coords = [this.points.coord];
            % vectorize elements
            this.elements(this.nEl) = Face;
            count = 1;
            for i = 1:length(this.elementsLevel)
                els = this.elementsLevel{i};
                this.elements(count:count+length(els)-1) = els;
                count = count+length(els);
            end
            
        end
        function initializeEl(this)
            n = length(this.elements);
            nMother = length(this.elementsMother);
            %The mother elements doesnt need to initialize
            for i = 1:nMother
                el = this.elements(i);
                el.initialize(this.points);
            end
            for i = nMother+1:n
                el = this.elements(i);
                el.initialize(this.points);
                %el.findAdjacents(this.heds,this.edges);
            end
        end
        function points = getPointsInLevel(this,nv)
            elements = this.elementsLevel{nv};
            ids = [];
            points = Point.empty;
            for i = 1:length(elements)
                el = elements(i);
                for j = 1:length(el.refinedPoints)
                    pi = el.refinedPoints(j);
                    if ~ismember(pi.id,ids)
                        pi.elements(end+1) = el;
                        points(end+1) = pi;
                        ids(end+1) = pi.id;
                    else
                        pi.elements(end+1) = el;
                    end
                end
                if strcmp(this.typeEl,'Const')
                    continue;
                end
                for j = 1:length(el.points)
                    pi = el.points(j);
                    if ~ismember(pi.id,ids)
                        pi.elements(end+1) = el;
                        points(end+1) = pi;
                        ids(end+1) = pi.id;
                    else
                        pi.elements(end+1) = el;
                    end
                end
            end
        end
        function initializePoints(this)
            elements = this.elementsLevel{end};
            for (element = elements)
                for (pt = element.geometry.geometryNodes())
                    pt.elements(end+1) = element;
                end
            end
        end
        function fig = plot(this,nv,color,alpha,q)
            if nargin == 2
                color = 'blue';
                alpha = 0.5;
            end
            if nargin == 3
                alpha = 0.5;
            end
            fig = figure;
            axis vis3d;
            axis equal
            axis off
            camproj('persp');
            T = get(gca,'tightinset');
            set(gcf,'color','white');
            lighting gouraud;
            hold on
            set(gca,'XColor', 'none','YColor','none','ZColor','none');
            set(gcf,'color','w');
            %fig.Position = [100 100 1024  768];
            view(34,41);
            elementsNv  = this.elementsLevel{nv};
            inc = reshape([elementsNv.pointsInc]',3,[])';
            pts = [this.points.coord]';
            if nargin == 5
                trisurf(inc,pts(:,1), pts(:,2), pts(:,3),q);
                shading interp
                colormap('jet');
            else
                trisurf(inc,pts(:,1), pts(:,2), pts(:,3),'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.95);
            end            
            this.Shadow(inc,pts(:,1), pts(:,2), pts(:,3));
            %set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);

            %             for i = 1:length(elementsNv)
            %                 el = elementsNv(i);
            %                 el.plot(color,alpha);
            %             end
            %get Centroide
            %pos = this.centroide + [8;0;0];
            %text(pos(1),pos(2),pos(3),strcat('NL =  ',num2str(nv)));
            %hold off;
        end
        function fig = plotByElement(this,nv, color,alpha)
            fig = figure;
            hold on
            axis equal;
            ax = gca;
            set(gca,'XColor', 'none','YColor','none','ZColor','none');
            set(gcf,'color','w');
            lighting gouraud;
            a1 = camlight('right');
            a2 = camlight('left');
            a3 = camlight('headlight');
            fig.Position = [100 100 1024  768];
            view(150,2);
            elementsNv  = this.elementsLevel{nv};
            pts = this.getPtsFromElements(elementsNv);
            for pt = pts
                pt.plot();
            end
            %count = 1;
            %exportgraphics(gca,strcat('element',int2str(count) ,'.png'),'Resolution',1000);
            for i = 1:length(elementsNv)
                el = elementsNv(i);
                el.plot(color,alpha);
            end
            % count = count + 1;
            %exportgraphics(gca,strcat('element',int2str(count) ,'.png'),'Resolution',1000);
            %             for i = 1:length(elementsNv)
            %                 el = elementsNv(i);
            %                 el.plotHed(this);
            %             end
            %             for edge = this.edges
            %                 edge.plot(this);
            %             end
            hold off
            
            
        end
        function pts = getPtsFromElements(this,elements)
            inc = [];
            for el = elements
                incEl = el.pointsInc;
                inc(end+1: end+3 ) = incEl;
            end
            inc = unique(inc);
            pts = this.points(inc);
        end
        
    end
    methods(Static)
        function Shadow(conec,x,y,z)
            nudge = 0.5;
            %l = light('Position',[1 4 5.0],'Style','infinite');
            l = light('Position',[1 1 5.0],'Style','infinite');
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
                D = Solid.matrixnormalize(D);
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
    end
end

