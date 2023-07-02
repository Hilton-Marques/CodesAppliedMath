classdef OctreeMethod < handle
    %This is used to generate a octree data structure.
    %   Properties:
    %       starPoints, starBinMatch, binCount, binBoundary, binDepth, 
    %       binParents, binCapacity, maxSize, minSize
    %   Methods:
    %       OctreeMethod, devide, devideBin, shrink
    %   
    %   THANK Sven for his work OcTree(2013) on MATLAB File Exchange.
    
    properties
        starPoints
        starBinMatch
        binCount
        binBoundary
        binDepth
    end
    properties(Hidden)
        binParents = zeros(0,1);
        binCapacity
        maxSize = inf;
        minSize = 1000*eps;
    end
    
    methods
        function obj = OctreeMethod(starCluster,starNumber) 
            pts = reshape([starCluster.position],3,starNumber)';
            obj.binBoundary = [min(pts,[],1) max(pts,[],1)];
            obj.starPoints = pts;
            obj.starBinMatch = ones(starNumber,1);
            obj.binDepth = 0;
            obj.binParents(1) = 0;
            obj.binCount = 1;
            obj.binCapacity = ceil(starNumber)/10;
            
            % Start dividing
            obj.preallocateSpace;
            obj.divide(1);
            obj.deallocateSpace;
        end
        
        function shrink(obj)
            % Shrink all bins to bound only the points they contain
            binChildren = arrayfun(@(i)find(obj.binParents==i),1:obj.binCount,'Un',0)';
            binIsLeaf = cellfun(@isempty, binChildren);
            for i = find(binIsLeaf(:))'
                binShrink_recurse(i, true)
            end
            
            function binShrink_recurse(binNo, isLeafBin)
                % Build a list of all points that fall within one of the
                % bins to be checked, and the list of which point falls in
                % which bin.
                oldBoundaryMin = obj.binBoundary(binNo,1:3);
                oldBoundaryMax = obj.binBoundary(binNo,4:6);
                if isLeafBin
                    % Shrink bin based on child POINTS
                    ptsMask = obj.starBinMatch==binNo;
                    if ~any(ptsMask)
                        % No points, shrink the bin to infinitely small
                        proposedBoundaries = [oldBoundaryMin oldBoundaryMin];
                    else
                        pts = obj.starPoints(ptsMask,:);
                        proposedBoundaries = [...
                            max([oldBoundaryMin; min(pts,[],1)]) ...
                            min([oldBoundaryMax; max(pts,[],1)])];
                    end
                else
                    % Shrink bin based on child BINS
                    childBoundaries = obj.binBoundary(binChildren{binNo},:);
                    proposedBoundaries = [min(childBoundaries(:,1:3),[],1) max(childBoundaries(:,4:6),[],1)];
                end
                
                if ~isequal(proposedBoundaries, [oldBoundaryMin oldBoundaryMax])
                    obj.binBoundary(binNo,:) = proposedBoundaries;
                    parentBin = obj.binParents(binNo);
                    if parentBin>0
                        binShrink_recurse(parentBin, false)
                    end
                end
            end
        end
    end
    
    methods(Access = protected)
        function preallocateSpace(obj)
            numPts = size(obj.starPoints,1);
            numBins = numPts;
            if isfinite(obj.binCapacity)
                numBins = ceil(2*numPts/obj.binCapacity);
            end
            obj.binDepth(numBins) = 0;
            obj.binParents(numBins) = 0;
            obj.binBoundary(numBins,1) = 0;
        end
        function deallocateSpace(obj)
            obj.binDepth(obj.binCount+1:end) = [];
            obj.binParents(obj.binCount+1:end) = [];
            obj.binBoundary(obj.binCount+1:end,:) = [];
        end
        
        function divide(obj, startingBins)
            % Loop over each bin we will consider for division
            for i = 1:length(startingBins)
                binNo = startingBins(i);
                
                % Prevent dividing beyond a minimum size                
                objBounds = obj.binBoundary(binNo,:);
                binEdgeSize = diff(objBounds([1:3;4:6]));
                minEdgeSize = min(binEdgeSize);
                maxEdgeSize = max(binEdgeSize);
                if minEdgeSize < obj.minSize
                    continue;
                end
                
                % There are two conditions under which we should divide
                % obj bin. 1: It's bigger than maxSize. 2: It contains
                % more points than binCapacity.
                oldCount = obj.binCount;
                if nnz(obj.starBinMatch==binNo) > obj.binCapacity
                    obj.divideBin(binNo);
                    obj.divide(oldCount+1:obj.binCount);
                    continue;
                end
                if maxEdgeSize > obj.maxSize
                    obj.divideBin(binNo);
                    obj.divide(oldCount+1:obj.binCount);
                    continue;
                end
            end
        end
        
        function divideBin(obj,binNo)
            % Gather the new points (a bit more efficient to copy once)
            binPtMask = obj.starBinMatch==binNo;
            objBinsPoints = obj.starPoints(binPtMask,:);
            
            % Get the old corner points and the new division point
            oldMin = obj.binBoundary(binNo,1:3);
            oldMax = obj.binBoundary(binNo,4:6);
            newDiv = mean([oldMin; oldMax], 1);
            
            % Build the new boundaries of our 8 subdivisions
            minMidMax = [oldMin newDiv oldMax];
            newBounds = minMidMax([...
                1 2 3 4 5 6;
                1 2 6 4 5 9;
                1 5 3 4 8 6;
                1 5 6 4 8 9;
                4 2 3 7 5 6;
                4 2 6 7 5 9;
                4 5 3 7 8 6;
                4 5 6 7 8 9]);
            
            % Determine to which of these 8 bins each current point belongs
            binMap = cat(3,[0 0 0],[0 0 1],[0 1 0],[0 1 1],...
                [1 0 0],[1 0 1],[1 1 0],[1 1 1]);
            gtMask = bsxfun(@gt, objBinsPoints, newDiv);
            [~,binAssignment] = max(all(bsxfun(@eq,gtMask,binMap),2),[],3);
            
            % Make the new bins and reassign old points to them
            newBinInds = obj.binCount+1:obj.binCount+8;
            obj.binBoundary(newBinInds,:) = newBounds;
            obj.binDepth(newBinInds) = obj.binDepth(binNo)+1;
            obj.binParents(newBinInds) = binNo;
            obj.starBinMatch(binPtMask) = newBinInds(binAssignment);
            obj.binCount = obj.binCount + 8;
        end
    end
end

