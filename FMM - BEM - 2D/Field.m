classdef Field < handle
    properties
        root
        nCells
        bodies
        cellHash = {};
        levels = [];
        maxLevel;
        k
        ax
        ng = 30;
    end
    methods(Static)
        function [side,centroid] = findConvexHull(posVec)
            boundary = [min(posVec,[],1) max(posVec,[],1)];
            centroid = mean([boundary(1:2); boundary(3:4)], 1);
            side = max((boundary(3) - boundary(1)),(boundary(4) - boundary(2)));
        end
        %% Get Gauss points
        function [x,w]= GaussQuad(N,a,b)
            N = N-1;
            N1 = N+1; N2=N+2;
            xu = linspace(-1,1,N1)';
            % Initial guess
            y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
            % Legendre-Gauss Vandermonde Matrix
            L=zeros(N1,N2);
            % Derivative of LGVM
            Lp=zeros(N1,N2);
            % Compute the zeros of the N+1 Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method
            y0=2;
            % Iterate until new points are uniformly within epsilon of old points
            while max(abs(y-y0))>eps
                L(:,1)=1;
                Lp(:,1)=0;
                L(:,2)=y;
                Lp(:,2)=1;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                y0=y;
                y=y0-L(:,N2)./Lp;
            end
            % Linear map from[-1,1] to [a,b]
            x=(a*(1-y)+b*(1+y))/2;
            % Compute the weights
            w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
        end
    end
    methods
        function this = Field(nodes,capacity,k)
            this.k = k;
            this.ax = zeros(length(nodes),1);
            this.bodies = nodes;
            pos = [nodes.position];
            posVec = [ [pos(1:2:end)]' [pos(2:2:end)]'];
            [side,centroid] = this.findConvexHull(posVec);
            side = 7.678125;
            this.root = Leaf(side,centroid,capacity,'',0,k);
            this.start();
            [this.cellHash,this.levels] = this.root.buildCellHash(this.cellHash,this.levels);
            this.maxLevel = size(this.levels,2);
            this.root.showAll();
            tic
            %% Upward
            for i = this.maxLevel:-1:3
                index = logical(this.levels(:,i));
                cells = this.cellHash(index,3);
                for j = 1:length(cells)
                    cell = cells{j};
                    if (~cell.isDivided)
                        %Eq. 3.19
                        this.calculateME(cell);
                    else
                        %Eq. 3.23
                        this.calculateMMT(cell)
                    end
                end
            end
            %% Doward Pass
            for i = 3:this.maxLevel
                index = logical(this.levels(:,i));
                cells = this.cellHash(index,3);
                for j = 1:length(cells)
                    cell = cells{j};
                    iteractionList = this.findInteractionList(cell.key);
                    for k = 1:length(iteractionList)
                        field = iteractionList(k);
                        %Eq. 3.25
                        this.calculateML(cell,field);
                    end
                    if (cell.isDivided())
                        %Eq. 3.30
                        this.calculateLLT(cell);
                    else
                        this.calculateNearInt(cell);
                        this.calculateFarInt(cell);
                    end
                end
            end
            this.ax
        end
        %% Create root
        function start(this)
            bodies_ = this.bodies;
            for i = 1:length(bodies_)
                this.root.insertBody(bodies_(i));
            end
        end
        %% Find Neighborhood
        function out = findVicinity(this,key)
            level = numel(key)+1;
            indexLevel = logical(this.levels(:,level));
            cells = this.cellHash(indexLevel,2:3);
            vicinityIndex = this.vicinityIndex(key);
            chaves = cells(:,1);
            chavesDec = base2dec(chaves,4);
            [index,~]=ismember(chavesDec,vicinityIndex);
            out = cells(index,2);
            out = [out{:}];
        end
        function out = findAdjacentCells(this,key)
            out = this.findVicinity(key);
            parent = key(1:end-1);
            candidatos = this.findVicinity(parent);
            for i = 1:length(candidatos)
                cell = candidatos(i);
                if (~cell.isDivided)
                    out(end+1) = cell;
                end
            end
        end
        function out = findInteractionList(this,key)
            parent = key(1:end-1);
            iteractionList = {};
            parents = this.findVicinity(parent);
            if (isempty(parents))
                out = [];
                return;
            end
            parents = parents([parents(:).isDivided]);
            keyVicinity = this.vicinityIndex(key);
            for i = 1:length(parents)
                cellP = parents(i);
                children = [cellP.children{:}];
                children = children([children(:).nBodies] > 0);
                for j = 1:length(children)
                    cell_ = children(j);
                    iteractionList(end+1,1:2) = {base2dec(cell_.key,4),...
                        cell_};
                end
            end
            if isempty(iteractionList)
                out = [];
                return;
            end
            [index,~] = ismember([iteractionList{:,1}],keyVicinity);
            out = {iteractionList{~index,2}};
            out = [out{:}];
        end
        %% Morton Code
        function out = interleave(~,bin)
            for i = 1:length(bin)
                be = bin(i);
                beBin = dec2bin(be);
                binOdd = strcat('0b',beBin(2:2:end));
                binEven = strcat('0b',beBin(1:2:end));
                out(i,1:2) = [str2num(binOdd),str2num(binEven)];
            end
        end
        function out = decode(~,bin)
            binStr = dec2bin(bin);
            if (mod(numel(binStr),2) == 1)
                binStr = strcat('0',binStr);
            end
            binOdd = strcat('0b',binStr(2:2:end));
            binEven = strcat('0b',binStr(1:2:end));
            out = [int16(str2num(binOdd)),int16(str2num(binEven))];
        end
        function out = encode(~,bin)
            x = bin(:,1);
            y = bin(:,2);
            xBin = dec2bin(x);
            yBin = dec2bin(y);
            nEl = max(numel(xBin(1,:)),numel(yBin(1,:)));
            xBin = dec2bin(x,nEl);
            yBin = dec2bin(y,nEl);
            out = char(zeros(length(x),2*nEl));
            out(:,2:2:end) = xBin;
            out(:,1:2:end) = yBin;
            out = strcat('0b',out);
            out = str2num(out);
        end
        function out = vicinityIndex(this,key)
            keyDec = base2dec(key,4);
            keyXY = this.decode(keyDec);
            vicinity = [(keyXY(1)-1),keyXY(1),(keyXY(1)+1),(keyXY(1)+1),...
                (keyXY(1)+1),keyXY(1),(keyXY(1)-1),(keyXY(1) -1);(keyXY(2)-1),...
                (keyXY(2)-1),(keyXY(2)-1),keyXY(2),(keyXY(2)+1),(keyXY(2)+1),...
                (keyXY(2)+1),keyXY(2)];
            out = this.encode(vicinity');
        end
        function parents = findParents(this,leaf,parents)
            if (leaf.level <= 2)
                return
            end
            parentKey = leaf.key(1:end-1);
            mother = this.cellHash{strcmp(this.cellHash(:,2),parentKey),3};
            parents{end+1} = mother;
            parents = this.findParents(mother,parents);
        end
        
        %% Math Calculus
        function calculateME(~,leaf)
            zc = complex(leaf.nodeCE(1),leaf.nodeCE(2));
            for i = 1:leaf.nBodies
                body = leaf.bodies(i);
                element = body.element;
                za = complex(element.position(1,1),element.position(2,1));
                zb = complex(element.position(3,1),element.position(4,1));
                tan = (zb - za)/norm(zb-za);
                n = -complex(0,1)*tan;
                wb = conj(tan);
                z1 = za - zc;
                z2 = zb - zc;
                zp1 = wb*z1;
                zp2 = wb*z2;
                if (element.q(1) == 1)
                    phi = element.q(2);
                    q = 0;
                else
                    q = element.q(2);
                    phi = 0;
                end
                leaf.MEF(1) = leaf.MEF(1) +  -(zp2 - zp1)*q;
                for j = 2:leaf.k+1
                    leaf.MEF(j) = leaf.MEF(j) + phi*n*(zp2 - zp1);
                    zp1 = zp1*z1/(j);
                    zp2 = zp2*z2/(j);
                end
                %this.showVector(leaf.nodeCE',[0.6,0.2,0.2],body.position');
            end
        end
        function calculateMMT(~,cell)
            children = [cell.children{:}];
            children = children([children(:).nBodies] > 0);
            zcP = complex(cell.nodeCE(1),cell.nodeCE(2));
            for i = 1:length(children)
                child = children(i);
                zc = complex(child.nodeCE(1),child.nodeCE(2));
                z0 = zc - zcP;
                k_ = child.k;
                zi = complex(1,0);
                for j = 1:k_+1
                    for m = j:k_+1
                        cell.MEG(m,:) = cell.MEG(m,:)+ zi*child.MEG(m-j+1,:);
                        cell.MEF(m,:) = cell.MEF(m,:) + zi*child.MEF(m-j+1,:);
                    end
                    zi = zi*z0/(j);
                end
            end
        end
        function calculateML(~,ceSource,ceField)
            p = ceField.k;
            zSource = complex(ceSource.nodeCE(1),ceSource.nodeCE(2));
            zField = complex(ceField.nodeCE(1),ceField.nodeCE(2));
            z0 = zSource - zField;
            zo = 1;
            ceSource.ML(1) = ceSource.ML(1) - log(z0)*ceField.MEF(1);
            for m = 2: 2*(p+1)
                zo = zo/z0;
                kmin = max(1,m - p);
                kmax = min(m,p+1);
                sgn = (-1)^(kmin-1);
                for j = kmin:kmax
                    ceSource.ML(j) = ceSource.ML(j) + ...
                        sgn*zo*ceField.MEF(m-j+1);
                    sgn = - sgn;
                end
                zo = zo*(m-1);
            end
            %this.showVector(ceSource',[0.6,0.2,0.2],ceField');
        end
        function calculateLLT(~,cell)
            zCp = complex(cell.nodeCE(1),cell.nodeCE(2));
            children = [cell.children{:}];
            children = children([children(:).nBodies] > 0);
            for i = 1:length(children)
                child = children(i);
                zCc = complex(child.nodeCE(1),child.nodeCE(2));
                z0 = (zCc - zCp);
                zi = 1;
                for j = 1:cell.k+1
                    for m = 1:(cell.k+1) - j + 1
                        child.ML(m) = child.ML(m) + zi*cell.ML(j+m-1,:);
                    end
                    zi = zi*z0/(j);
                end
            end
        end
        function calculateNearInt(this,source)
            bodiesSource = source.bodies;
            for i = 1:length(bodiesSource)
                bodyS = bodiesSource(i);
                elementS = bodyS.element;
                for j = 1:length(bodiesSource)
                    bodyF = bodiesSource(j);
                    elementF = bodyF.element;
                    if (elementF.q(1,1) == 1)
                        if (i == j)
                            this.ax(elementS.id) = this.ax(elementS.id) + ...
                                0.5*elementF.q(2,1);
                        else
                            this.ax(elementS.id) = this.ax(elementS.id) - ...
                                this.calculateHHd(bodyS.position,elementF);
                        end
                    else
                        if (i == j)
                            this.ax(elementS.id) = this.ax(elementS.id) + elementF.q(2,1)*...
                                ((elementF.len)/(pi*2))*(log(elementF.len/2)-1);
                        else
                            this.ax(elementS.id) = this.ax(elementS.id) - ...
                                this.calculateGGq(bodyS.position,elementF);
                        end
                    end
                end
                fields = this.findAdjacentCells(source.key);
                for j = 1:length(fields)
                    field = fields(j);
                    bodiesF = field.bodies;
                    for l = 1:length(bodiesF)
                        bodyF = bodiesF(l);
                        elementF = bodyF.element;
                        if (elementF.q(1,1) == 1)
                            this.ax(elementS.id) = this.ax(elementS.id) - ...
                                this.calculateHHd(bodyS.position,elementF);
                        else
                            this.ax(elementS.id) = this.ax(elementS.id) - ...
                                this.calculateGGq(bodyS.position,elementF);
                        end
                    end
                end
            end
        end
        function out = calculateHHd(this,xo,elementF)
            out = 0;
            [ksi,w] = this.GaussQuad(this.ng,0,1);
            for i = 1:this.ng
                xksi = elementF.position(1) + (elementF.position(3) ...
                    - elementF.position(1))*ksi(i);
                yksi = elementF.position(2) + (elementF.position(4) ...
                    - elementF.position(2))*ksi(i);
                r = [xksi - xo(1); yksi - xo(2)];
                n = [elementF.position(4) - elementF.position(2); ...
                    -(elementF.position(3) - elementF.position(1))];
                n = n/vecnorm(n);
                out = out + (elementF.q(2,1)*(elementF.len/(2*pi*norm(r)^2))* ...
                    dot(r,n))*w(i);
            end
        end
        function out = calculateGGq(this,xo,elementF)
            out = 0;
            [ksi,w] = this.GaussQuad(this.ng,0,1);
            for i = 1:this.ng
                xksi = elementF.position(1) + (elementF.position(3) ...
                    - elementF.position(1))*ksi(i);
                yksi = elementF.position(2) + (elementF.position(4) ...
                    - elementF.position(2))*ksi(i);
                r = [xksi - xo(1); yksi - xo(2)];
                out = out + (elementF.len/(2*pi))*log(vecnorm(r))*w(i);
            end
        end
        function calculateFarInt(this,cell)
            bodies = cell.bodies;
            fact = 1;
            for i = 2:this.k+1
                fact = fact/(i-1);
                cell.ML(i) = cell.ML(i)*fact;
            end
            for i = 1:length(bodies)
                body = bodies(i);
                zp = cell.ML(end);
                %z0 = zsource - zCentroid 
                z0 = complex(body.position(1) - cell.nodeCE(1),...
                    body.position(2) - cell.nodeCE(2));
                for j = this.k:-1:1
                    zp = zp*z0 + cell.ML(j);
                end
                zp = zp/(2*pi);
                this.ax(body.element.id) = this.ax(body.element.id) + real(zp);
            end
        end
        %% Plot
        function show(~,leafs,color)
            if (isempty(leafs))
                return
            end
            if (nargin < 3)
                color = [1 0 0];
            end
            for i = 1:length(leafs)
                leafs(i).showCurrent(color);
            end
        end
        function out = showVector(this,xf,c,xi)
            if (nargin < 2)
                c = [1 0 1];
            end
            if (nargin < 3)
                xi = zeros(2,1);
            end
            X = [xi(1,1), xf(1,1)];
            Y = [xi(2,1), xf(2,1)];
            h1 = line( X, Y, 'Color', c,'Linewidth',2);
            h2 = this.triangle(xf,c,xi);
            h = [h1,h2];
            pause(1);
            delete(h);
        end
        function out = triangle(~,xf,c,xi)
            x = xf - xi;
            L = vecnorm(x);
            factorReduc = 1/5;
            h = factorReduc*0.25;
            b = factorReduc*0.25;
            p_x = [-x(2,1); x(1,1)];
            p_x = p_x/norm(p_x);
            p1 = xf + b/2*p_x;
            p2 = xf - b/2*p_x;
            p3 = xf + h*(x/norm(x));
            X = [p1(1,1), p2(1,1) , p3(1,1)];
            Y = [p1(2,1), p2(2,1) , p3(2,1)];
            out = fill(X, Y, c);
        end
        %% Find square whose bodies are inside
        
    end
end