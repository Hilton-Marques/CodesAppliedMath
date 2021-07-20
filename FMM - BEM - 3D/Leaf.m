classdef Leaf < handle
    properties
        % OctTree variables
        s;  %side
        nodeCE = []; %leaf's centroid
        nodes = Node.empty; % nodes beloging
        capacity; % max capacity
        isDivided = false;
        nNodes = 0;
        children = cell(8,1);
        nChildren = [];
        level;
        key;
        % FMM variables
        n;   % max expansion
        MEG; % moments G matrix
        MEH; % moments H matrix 
        MLG; % moments local G
        MLH; % moments local H
    end
    methods
        function this = Leaf(s,cellCentroid,capacity,key,level,n)
            this.s = s;
            this.nodeCE = cellCentroid;
            this.capacity = capacity;
            this.key = key;
            this.level = level;
            this.n = n;
            % Total of terms 
            a1 = 1;
            an = 2*(n)+1;
            totN = (a1 + an)*(n+1)/2;
            this.MEG = zeros(totN,1);
            this.MEH = zeros(totN,1);
            this.MLG = zeros(totN,1);
            this.MLH = zeros(totN,1);
        end
        function insertNode(this,node)
            if this.nNodes < this.capacity
                this.update(node);
            else
                if(~this.isDivided)
                    this.subdivide();
                end
                this.update(node);
                newOctant = this.getOctant(node);
                this.children{newOctant}.insertNode(node);
            end
        end
        function update(this,node)
            this.nNodes = this.nNodes + 1;
            this.nodes(this.nNodes) = node;
        end
        function subdivide(this)
            points = [this.nodes.pos];
            position = [ [points(1:3:end)]' [points(2:3:end)]' [points(3:3:end)]'];
            gtMidMask = position >= this.nodeCE;
            octant = gtMidMask(:,1) + 2*gtMidMask(:,2) + 4*gtMidMask(:,3) + 1;
            this.setChild(octant);
            this.isDivided = true;
        end
        function setChild(this,octant)
            s = this.s;
            oldCE = this.nodeCE;
            capacity = this.capacity;
            idxOctant = zeros(length(octant),8);
            nNodesbyChild = zeros(8,1);
            for i = 1:8
                idxOctant(:,i) = octant == i;
                nNodesbyChild(i) = sum(idxOctant(:,i));
            end
            newLevel = this.level + 1;
            this.children{1} = Leaf(s/2,oldCE + [-s/4,-s/4,-s/4],capacity,...
                strcat(this.key,'000'),newLevel,this.n);
            this.children{2} = Leaf(s/2,oldCE +  [s/4,-s/4,-s/4],capacity,...
                strcat(this.key,'001'),newLevel,this.n);
            this.children{3} = Leaf(s/2,oldCE +  [-s/4,s/4,-s/4],capacity,...
                strcat(this.key,'010'),newLevel,this.n);
            this.children{4} = Leaf(s/2,oldCE +  [s/4,s/4,-s/4] ,capacity,...
                strcat(this.key,'011'),newLevel,this.n);
            this.children{5} = Leaf(s/2,oldCE + [-s/4,-s/4,s/4],capacity,...
                strcat(this.key,'100'),newLevel,this.n);
            this.children{6} = Leaf(s/2,oldCE +  [s/4,-s/4,s/4],capacity,...
                strcat(this.key,'101'),newLevel,this.n);
            this.children{7} = Leaf(s/2,oldCE +  [-s/4,s/4,s/4],capacity,...
                strcat(this.key,'110'),newLevel,this.n);
            this.children{8} = Leaf(s/2,oldCE +  [s/4,s/4,s/4] ,capacity,...
                strcat(this.key,'111'),newLevel,this.n);
            for i = 1:8
                newNNodes = nNodesbyChild(i); % new Nodes beloging to chilg
                if (newNNodes > 0)
                    this.children{i}.nodes(1:newNNodes,1) = ...
                        [this.nodes(logical(idxOctant(:,i)))];
                    this.children{i}.nNodes = newNNodes;
                    this.nChildren(end+1) = 1;
                end
            end
        end
        function output = getOctant(this,node)
            pos = node.pos;
            gtMidMask = pos >= this.nodeCE;
            output = gtMidMask(:,1) + 2*gtMidMask(:,2) + 4*gtMidMask(:,3) + 1;
        end
        function [cellHash,levels] = buildCellHash(this,cellHash,levels)
            cellkey = this.key;
            cellHash(end+1,1:2) = {cellkey,this};
            levels(end+1,(numel(cellkey)-2)/3 + 1) = 1;
            if (~this.isDivided)
                return;
            end
            for i = 1:8
                child = this.children{i};
                if (child.nNodes > 0)
                    [cellHash,levels] = child.buildCellHash(cellHash,levels);
                end
            end
        end
        function showAll(this)
            left_corner_x = this.nodeCE(1,1) - this.s*0.5;
            left_corner_y = this.nodeCE(1,2) - this.s*0.5;
            rectangle('Position',[left_corner_x,left_corner_y,this.s,this.s]);
            if (this.isDivided)
                for i = 1:8
                    this.children{i,1}.showAll();
                end
            end
        end
    end
end