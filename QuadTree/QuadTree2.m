classdef QuadTree2 < handle
    properties
        points = [];
        idPointsToCells = [];
        boundary = [];
        cell_level = [];
        cell_parent = [];
        cell_count = [];
        capacity = [];
    end
    methods
        function this = QuadTree2(pts,capacity)
            numPts = size(pts,1);
            boundary = [min(pts,[],1) max(pts,[],1)];
            w = 1.1*(boundary(3) - boundary(1));
            h = 1.1*(boundary(4) - boundary(2));
            cellCentroid = mean([boundary(1:2); boundary(3:4)], 1);
            this.boundary = [cellCentroid - [w/2,h/2], cellCentroid + [w/2,h/2]];
            this.points = pts;
            this.cell_count = 1;
            this.cell_level = 0;
            this.cell_parent = 0;
            this.idPointsToCells = ones(numPts,1);
            this.capacity = capacity;
            this.divide(1);
        end
        function divide(this,startingLevel)
            %Loop for each cell we will consider for division
            for i= 1:length(startingLevel)
                cellId = startingLevel(i);
                oldCellsNumber = this.cell_count;
                %ask: how many points are inside a given cell
                if nnz(this.idPointsToCells == cellId) > this.capacity
                    this.divideCell(cellId);
                    this.divide(oldCellsNumber+1:this.cell_count);
                    continue;
                end
            end
        end
        function divideCell(this,cellId)
            %find: wich points belongs to a given cell
            points_id = this.idPointsToCells == cellId; %True or false
            cellPoints = this.points(points_id,:);
            
            %Get the dimensions of a cell
            cellBoundMin = this.boundary(cellId,1:2);
            cellBoundMax = this.boundary(cellId,3:4);
            cellCentroid = mean([cellBoundMin; cellBoundMax], 1);
            cellRec = [cellBoundMin cellCentroid cellBoundMax];
            childrenBounds = cellRec([...
                1 2 3 4;...
                3 2 5 4;...
                1 4 3 6;...
                3 4 5 6]);
            %Map southwest,southeast,northwest,northeast to logical
            cellMap = [ [0 0]; [1 0]; [0 1]; [1 1]];%1->Sw,2->SE,3->NW,4->NE
            gtMidMask = cellPoints > cellCentroid;
            %compare both and get index
            [~,idPointsToChildren] = ismember(gtMidMask,cellMap,'rows');
            %Verify if all the children are populated
            newCellsId = [ismember([1;2;3;4],idPointsToChildren)];
            newCellsName = this.cell_count + newCellsId + [0;newCellsId(1:end-1)]+...
                [0;0;newCellsId(1:2)] + [0;0;0;newCellsId(1)];
            this.boundary(newCellsName(newCellsId),:) = childrenBounds(newCellsId,:);
            this.cell_level(newCellsName(newCellsId)) = this.cell_level(cellId) + 1;
            this.cell_parent(newCellsName(newCellsId)) = cellId;
            this.idPointsToCells(points_id) = newCellsName(idPointsToChildren);
            this.cell_count = this.cell_count + sum(newCellsId);
        end
        function show(this)
            hold on;
            for i = 1:this.cell_count
                boundary = this.boundary(i,:);
                w = boundary(3) - boundary(1);
                h = boundary(4) - boundary(2);
                left_corner = boundary([1 2]);
                rectangle('Position',[left_corner(1),left_corner(2),w,h]);
            end
        end
    end
end