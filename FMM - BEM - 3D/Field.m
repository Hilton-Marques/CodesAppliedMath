classdef Field < handle
    properties
        root
        nCells
        nodes
        cellHash = {};
        levels = [];
        maxLevel;
        n
        solver
    end
    methods
        function this = Field(nodes,capacity,n,ng)
            this.n = n;
            this.nodes = nodes;
            this.solver = Solver(ng,n,length(nodes));
            posVec = [nodes.pos];
            posVec = [ [posVec(1:3:end)]',[posVec(2:3:end)]',[posVec(3:3:end)]'];
            [side,centroid] = this.findBoundaryCube(posVec);
            this.root = Leaf(side,centroid,capacity,'0b',0,n);
            this.start();  % Start Octree
            [this.cellHash,this.levels] = this.root.buildCellHash(this.cellHash,this.levels);
            this.maxLevel = size(this.levels,2);
            %% Upward
            for i = this.maxLevel:-1:3
                index = logical(this.levels(:,i));
                cells = this.cellHash(index,2);
                for j = 1:length(cells)
                    cell = cells{j};
                    if (~cell.isDivided)
                        % Calculate multipole expansion Eq. 3.53 e 3.54
                        this.solver.calculateME(cell);
                    else
                        % Eq. 3.55
                        this.solver.calculateMMT(cell)
                    end
                end
            end
            %% Doward Pass
            for i = 3:this.maxLevel
                index = logical(this.levels(:,i));
                cells = this.cellHash(index,2);
                for j = 1:length(cells)
                    cell = cells{j};
                    iteractionList = this.findInteractionList(cell.key);
                    for k = 1:length(iteractionList)
                        field = iteractionList(k);
                        %Eq. 3.57
                        this.solver.calculateM2L(cell,field);
                    end
                    if (cell.isDivided())
                        %Eq. 3.58
                        this.solver.calculateLLT(cell);
                    else
                        % Eq. 3.56
                        this.solver.calculateFarInt(cell);
                        fields = this.findAdjacentCells(cell.key);
                        this.solver.calculateNearInt(cell,fields);
                    end
                end
            end
            b = real(this.solver.Ax);
            fileID = fopen('u.txt','w');
            fprintf(fileID,'%12.9f\n',b);
            fclose(fileID);
        end
        %% Create root
        function start(this)
            nodes_ = this.nodes;
            for i = 1:length(nodes_)
                this.root.insertNode(nodes_(i));
            end
        end
        %% OctTree Topology
        function out = findAdjacentCells(this,key)
            out = this.findVicinity(key);
            % Besides the vicinity we should look for the parents that
            % share a vertice with the current cell
            parent = key(1:end-3);
            candidatos = this.findVicinity(parent);
            for i = 1:length(candidatos)
                cell = candidatos(i);
                if (~cell.isDivided)
                    out(end+1) = cell;
                end
            end
        end
        function out = findInteractionList(this,key)
            parent = key(1:end-3);
            iteractionList = {};
            parents = this.findVicinity(parent);
            if (isempty(parents))
                out = [];
                return;
            end
            parents = parents([parents(:).isDivided]);
            vicinityKeys = this.findVicinityKeys(key);
            for i = 1:length(parents)
                cellP = parents(i);
                children = [cellP.children{:}];
                children = children([children(:).nNodes] > 0);
                for j = 1:length(children)
                    cell_ = children(j);
                    iteractionList(end+1,1:2) = {double(str2num(cell_.key)),cell_};
                end
            end
            if isempty(iteractionList)
                out = [];
                return;
            end
            [index,~] = ismember([iteractionList{:,1}]',vicinityKeys);
            out = {iteractionList{~index,2}};
            out = [out{:}];
        end
        %% Find Neighborhood
        function out = findVicinity(this,key)
            level = (numel(key)-2)/3 + 1;
            indexLevel = logical(this.levels(:,level));
            cells = this.cellHash(indexLevel,:);
            nCells = length(cells);
            keys = zeros(nCells,1);
            for i = 1:nCells
                keyBin = cells{i,1};
                keys(i) = str2num(keyBin);
            end
            vicinityKeys = this.findVicinityKeys(key);
            [index,~] = ismember(keys,vicinityKeys);
            out = cells(index,2);
            out = [out{:}];
        end
        function out = findVicinityKeys(this,key)
            keyXYZ = this.decode(key);
            % 26 neighboorhood
            vicinityBase = [(keyXYZ(1)-1),keyXYZ(1),(keyXYZ(1)+1),(keyXYZ(1)+1),...
                (keyXYZ(1)+1),keyXYZ(1),(keyXYZ(1)-1),(keyXYZ(1)-1),(keyXYZ(1));...
                (keyXYZ(2)-1),(keyXYZ(2)-1),(keyXYZ(2)-1),keyXYZ(2),...
                (keyXYZ(2)+1),(keyXYZ(2)+1),(keyXYZ(2)+1),keyXYZ(2),keyXYZ(2);
                (keyXYZ(3)-1), (keyXYZ(3)-1), (keyXYZ(3)-1), (keyXYZ(3)-1),(keyXYZ(3)-1), ...
                (keyXYZ(3)-1),(keyXYZ(3)-1),(keyXYZ(3)-1),(keyXYZ(3)-1)];
            vicinityMiddle = [(keyXYZ(1)-1),keyXYZ(1),(keyXYZ(1)+1),(keyXYZ(1)+1),...
                (keyXYZ(1)+1),keyXYZ(1),(keyXYZ(1)-1),(keyXYZ(1)-1);...
                (keyXYZ(2)-1),(keyXYZ(2)-1),(keyXYZ(2)-1),keyXYZ(2),...
                (keyXYZ(2)+1),(keyXYZ(2)+1),(keyXYZ(2)+1),keyXYZ(2);
                (keyXYZ(3)), (keyXYZ(3)), (keyXYZ(3)), (keyXYZ(3)),(keyXYZ(3)), ...
                (keyXYZ(3)),(keyXYZ(3)),(keyXYZ(3))];
            vicinitySuper = [(keyXYZ(1)-1),keyXYZ(1),(keyXYZ(1)+1),(keyXYZ(1)+1),...
                (keyXYZ(1)+1),keyXYZ(1),(keyXYZ(1)-1),(keyXYZ(1)-1),(keyXYZ(1));...
                (keyXYZ(2)-1),(keyXYZ(2)-1),(keyXYZ(2)-1),keyXYZ(2),...
                (keyXYZ(2)+1),(keyXYZ(2)+1),(keyXYZ(2)+1),keyXYZ(2),keyXYZ(2);
                (keyXYZ(3)+1), (keyXYZ(3)+1), (keyXYZ(3)+1), (keyXYZ(3)+1),(keyXYZ(3)+1), ...
                (keyXYZ(3)+1),(keyXYZ(3)+1),(keyXYZ(3)+1),(keyXYZ(3)+1)];
            vicinity = [vicinityBase,vicinityMiddle,vicinitySuper];
            out = this.encode(vicinity);
        end
        function out = decode(~,key)
            binx = strcat('0b',key(:,5:3:end));
            biny = strcat('0b',key(:,4:3:end));
            binz = strcat('0b',key(:,3:3:end));
            out = [int16(str2num(binx)),int16(str2num(biny)),int16(str2num(binz))];
        end
        function out = encode(~,bin)
            x = bin(1,:);
            y = bin(2,:);
            z = bin(3,:);
            xBin = dec2bin(x);
            yBin = dec2bin(y);
            zBin = dec2bin(z);
            nEl = max([numel(xBin(1,:)),numel(yBin(1,:)),numel(zBin(1,:))]);
            xBin = dec2bin(x,nEl);
            yBin = dec2bin(y,nEl);
            zBin = dec2bin(z,nEl);
            out = char(zeros(26,3*nEl));
            out(:,3:3:end) = xBin;
            out(:,2:3:end) = yBin;
            out(:,1:3:end) = zBin;
            out = str2num(strcat('0b',out));
        end
    end
    methods(Static)
        function [side,centroid] = findBoundaryCube(posVec)
            boundary = [min(posVec,[],1) max(posVec,[],1)];
            centroid = mean([boundary(1:3); boundary(4:6)], 1);
            side = max((boundary(3) - boundary(1)),(boundary(4) - boundary(2)));
        end
    end
end