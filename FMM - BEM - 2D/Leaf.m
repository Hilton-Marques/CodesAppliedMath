classdef Leaf < handle
    properties
        s;
        nodeCE = [];
        bodies = Node.empty;
        capacity = 0;
        isDivided = false;
        nBodies = 0;
        children = cell(4,1);
        nChildren = [];
        level;
        key;
        MEG;
        MEF;
        k;
        ML;
    end
    methods
        function this = Leaf(s,cellCentroid,capacity,key,level,k)
            %             this.bodies = [this.bodies;body];
            %             this.level = level;
            this.s = s;
            this.nodeCE = cellCentroid;
            this.capacity = capacity;
            %PreAllocate bodies
            this.key = key;
            this.level = level;
            this. k = k;
            this.MEF = zeros(k+1,1);
            this.MEG = zeros(k+1,1);
            this.ML = zeros(k+1,1);
        end
        function insertBody(this,body)
            if this.nBodies < this.capacity
                this.update(body);
            else
                if(~this.isDivided)
                    this.subdivide();
                end
                this.update(body);
                newQuadrant = this.getQuadrant(body);
                this.children{newQuadrant,1}.insertBody(body);
            end
        end
        function subdivide(this)
            points = [this.bodies.position];
            position = [ [points(1:2:end)]' [points(2:2:end)]'];
            gtMidMask = position >= this.nodeCE;
            quadrant = gtMidMask(:,1) + 2*gtMidMask(:,2) + 1;
            this.setChild(quadrant);
            this.isDivided = true;
        end
        function setChild(this,quadrant)
            s = this.s;
            oldCE = this.nodeCE;
            capacity = this.capacity;
            idx1quadrant = quadrant == 1;
            idx2quadrant = quadrant == 2;
            idx3quadrant = quadrant == 3;
            idx4quadrant = quadrant == 4;
            nBody1 = sum(idx1quadrant);
            nBody2 = sum(idx2quadrant);
            nBody3 = sum(idx3quadrant);
            nBody4 = sum(idx4quadrant);
            newLevel = this.level + 1;
            this.children{1,1} = Leaf(s/2,oldCE + [-s/4,-s/4],capacity,...
                strcat(this.key,'0'),newLevel,this.k);
            this.children{2,1} = Leaf(s/2,oldCE +  [s/4,-s/4],capacity,...
                strcat(this.key,'1'),newLevel,this.k);
            this.children{3,1} = Leaf(s/2,oldCE +  [-s/4,s/4],capacity,...
                strcat(this.key,'2'),newLevel,this.k);
            this.children{4,1} = Leaf(s/2,oldCE +  [s/4,s/4] ,capacity,...
                strcat(this.key,'3'),newLevel,this.k);
            if(nBody1 > 0)
                this.children{1,1}.bodies(1:nBody1,1) = [this.bodies(idx1quadrant)];
                this.children{1,1}.updateTotal(nBody1);
                this.nChildren(end+1) = 1;
            end
            if (nBody2 > 0)
                this.children{2,1}.bodies(1:nBody2,1) = [this.bodies(idx2quadrant)];
                this.children{2,1}.updateTotal(nBody2);
                this.nChildren(end+1) = 2;
            end
            if (nBody3 > 0)
                
                this.children{3,1}.bodies(1:nBody3,1) = [this.bodies(idx3quadrant)];
                this.children{3,1}.updateTotal(nBody3);
                this.nChildren(end+1) = 3;
            end
            if (nBody4 > 0)
                
                this.children{4,1}.bodies(1:nBody4,1) = [this.bodies(idx4quadrant)];
                this.children{4,1}.updateTotal(nBody4);
                this.nChildren(end+1) = 4;
            end
        end
        function output = getQuadrant(this,body)
            position = body.position';
            gtMidMask = position >= this.nodeCE;
            output = gtMidMask(:,1) + 2*gtMidMask(:,2) + 1;
        end
        function update(this,body)
            this.nBodies = this.nBodies + 1;
            this.bodies(this.nBodies,1) = body;
        end
        function updateTotal(this,nBodies)
            this.nBodies = nBodies;
        end
        function showAll(this)
            left_corner_x = this.nodeCE(1,1) - this.s*0.5;
            left_corner_y = this.nodeCE(1,2) - this.s*0.5;
            rectangle('Position',[left_corner_x,left_corner_y,this.s,this.s]);
            if (this.isDivided)
                for i = 1:4
                    this.children{i,1}.showAll();
                end
            end
        end
        function [h1,h2] = showCurrent(this,color,body)
            if nargin > 2
                plot(body.position(1,1),body.position(1,2),'o','Color','y');
            end
            left_corner_x = this.nodeCE(1,1) - this.s*0.5;
            left_corner_y = this.nodeCE(1,2) - this.s*0.5;
            h1 = rectangle('Position',[left_corner_x,left_corner_y,this.s,this.s],...
                'EdgeColor',color);
            h2 = plot(this.leafMassXCM(1,1),this.leafMassXCM(1,2),'o','Color','red');
        end
        function out = findLeaf(this,body)
            if (~this.isDivided)
                out = this;
                return
            end
            getQuadrant = this.getQuadrant(body);
            out = this.children{getQuadrant,1}.findLeaf(body);
        end
        function out = contains(this,pos)
            nodeCE = this.nodeCE;
            side = this.s;
            min = nodeCE - [s/2,s/2];
            max = nodeCE + [s/2,s/2];
            out = false;
            if ((pos(1,1) >= min(1,1) && pos(1,1) <= max(1,1)) && ...
                    (pos(2,1) >= min(2,1) && pos(2,1) <= max(2,1)))
                out = true;
            end
        end
        function out = vicinities(this,vicinityNodes,level)
            if (this.level ~= level)
                for i = 1:4
                    this.children{i,1}.vicinities(vicinityNodes,level);
                end
            elseif (ismember(this.nodeCE,vicinityNodes))
                out(end+1) = this;
            end
        end
        function [cellHash,levels] = buildCellHash(this,cellHash,levels)
            cellkey = this.key;
            cellHash(end+1,1:5) = {this.level,cellkey,this,~this.isDivided,...
                numel(cellkey) == 2};
            levels(end+1,numel(cellkey)+1) = 1;
            if (~this.isDivided)
                return;
            end
            for i = 1:4
                child = this.children{i,1};
                if (child.nBodies > 0)
                    [cellHash,levels] = child.buildCellHash(cellHash,levels);
                end
            end
        end
    end
end