classdef Field < handle
    properties
        root
        nCells
        bodies
        cellHash = {};
    end
    methods
        function this = Field(bodies,capacity)
            this.bodies = bodies;
            pos = [bodies.position];
            posVec = [ [pos(1:2:end)]' [pos(2:2:end)]'];
            [side,centroid] = this.findConvexHull(posVec);
            this.root = Leaf(side,centroid,capacity,'',0);
            this.start();
            this.cellHash = this.root.buildCellHash(this.cellHash);
            hold on
            this.root.showAll();
            %             for i = 1:length(bodies)
            %                 body = bodies(i);
            %                 leaf = this.root.findLeaf(body);
            %                 leaf.showCurrent([0 1 0]);
            %                 vicinity = this.findVicinity(leaf.key);
            %                 iteractionList = this.findInteractionList(leaf.key);
            %                 %this.show(vicinity);
            %                 %this.show(iteractionList,[0 0 1]);
            %             end
            %Upward Pass
            leafs = this.cellHash([this.cellHash{:,4}],3);
            for i = 1:length(leafs)
                leaf = leafs{i};
                parents = {};
                parents = this.findParents(leaf,parents);
                leaf.ME = this.calculateME(leaf);
                if (~isempty(parents))
                    for j = 1:length(parents)
                        parents{j}.ME = parents{j}.ME + leaf.ME;
                        this.showVector(parents{j}.nodeCE',[0.6,0.2,0.2],...
                            leaf.nodeCE');
                    end
                end
            end
            %Doward Pass
            cellsLevel2 = this.cellHash([this.cellHash{:,5}],3);
            for i = 1:length(cellsLevel2)
                cell = cellsLevel2{i};
                this.startDownard(cell);
            end
        end
        function [side,centroid] = findConvexHull(~,posVec)
            boundary = [min(posVec,[],1) max(posVec,[],1)];
            centroid = mean([boundary(1:2); boundary(3:4)], 1);
            side = max((boundary(3) - boundary(1)),(boundary(4) - boundary(2)));
        end
        function start(this)
            bodies_ = this.bodies;
            for i = 1:length(bodies_)
                this.root.insertBody(bodies_(i));
            end
        end
        function buildCellHash(this)
            root_ = this.root;
            cellHash_ = this.cellHash;
            count = 1;
            for i=1:this.maxLevel
                out = root_.getCellS(i);
                c = length(out);
                cellHash(count:c) = out;
                count = c + 1;
            end
        end
        function out = findVicinity(this,key)
            C = this.cellHash;
            Level = numel(key);
            vicinityIndex = this.vicinityIndex(key);
            chaves = [C{:,1}];
            COnLevel = C(chaves == Level,2:3);
            chavesDec = base2dec({COnLevel{:,1}},4);
            [index,~]=ismember(chavesDec,vicinityIndex);
            out = COnLevel(index,2);
            out = [out{:}];
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
        function startDownard(this,cell)
            iteractionList = this.findInteractionList(cell.key);
            for i = 1:length(iteractionList)
                field = iteractionList(i);
                cell.ML = cell.ML + this.calculateML(cell.nodeCE,field.nodeCE);
            end
            if (cell.isDivided())
                children = [cell.children{:}];
                children = children([cell(:).nBodies] > 0);
                for i = 1:length(children)
                    children(i).ML = children(i).ML + this.calculateML(children(i).nodeCE,...
                        cell.nodeCE);
                    this.startDownard(children(i));
                end
            else
                return;
            end
        end
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
        function parents = findParents(this,leaf,parents)
            if (leaf.level <= 2)
                return
            end
            parentKey = leaf.key(1:end-1);
            mother = this.cellHash{strcmp(this.cellHash(:,2),parentKey),3};
            parents{end+1} = mother;
            parents = this.findParents(mother,parents);
        end
        function out = calculateME(this,leaf)
            out = 0;
            for i = 1:leaf.nBodies
                body = leaf.bodies(i);
                r = leaf.nodeCE - body.position;
                out = out + dot(r,r);
                this.showVector(leaf.nodeCE',[0.6,0.2,0.2],body.position');
            end
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
        function out = calculateML(this,ceSource,ceField)
            r = (ceSource - ceField);
            out = dot(r,r);
            this.showVector(ceSource',[0.6,0.2,0.2],ceField');
        end
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
    end
end