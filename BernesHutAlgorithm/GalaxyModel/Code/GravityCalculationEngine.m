classdef GravityCalculationEngine
    %This is the part where gravitational calculation is executed.
    % Propertires:
    %   G, starCluster, starNumber, massStarArray, positionVectArray
    % Methods:
    %   GravityCalculationEngine, calculate, simpleCalculate, 
    %   octreeCalculate, isFar, binPosVect, addBinMass, binCheck
    %   - - - - - - - - - - - -
    %   Author: Haihong
    %           Valentine's Day, 2015
    %   inspired by Salman Khan's videos on astronomy at Khan Academy.
    
    properties(Hidden)
        G
        starCluster
        starNumber
        massStarArray
        positionVectArray
        octreeThreshold = 20;
        farThreshold = 0.2;
    end
    
    methods
        function obj = GravityCalculationEngine(CONST_G,starCluster,starNumber)
           obj.G = CONST_G;
           obj.starCluster = starCluster;
           obj.starNumber = starNumber;
           obj.massStarArray = [obj.starCluster.mass]';
           obj.positionVectArray = reshape([obj.starCluster.position],3,obj.starNumber)';
        end
        
        function starClusterNew = calculate(obj)
            if obj.starNumber < obj.octreeThreshold % Do not use octree. This is an adjustable threshold value.
                obj.simpleCalculate;
                starClusterNew = obj.starCluster;
            else % use octree
                obj.octreeCalculate;
                starClusterNew = obj.starCluster;
            end
        end
    end
    
    methods(Access = protected)
        function simpleCalculate(obj)
            for s = 1:obj.starNumber
                massStarArrayTemp = obj.massStarArray;
                massStarArrayTemp(s) = []; % Delete the s-th element, size changed
                posDiffVectArray = bsxfun(@minus,obj.positionVectArray,obj.starCluster(s).position);% starNumber rows, 3 columns
                posDiffVectArray(s,:) = []; % Delete the s-th row, which is [0,0,0]. Size changed
                distCubeArray = (sum(posDiffVectArray.^2,2)).^1.5; % sum(A,2) adds A's elements horizontally
                
                % Vectorized Universal Law of Gravitation.
                % Note that massArrayTemp and distCubeArray are vertical.
                sjAccVectNumerator = obj.G * bsxfun(@times,massStarArrayTemp,posDiffVectArray);
                sjAccVect = bsxfun(@rdivide,sjAccVectNumerator,distCubeArray);
                obj.starCluster(s).changeAcceleration(sum(sjAccVect,1)); % sum(A,1) adds A's elements 'vertically'
            end
        end
        
        function octreeCalculate(obj)
            octree = OctreeMethod(obj.starCluster,obj.starNumber);
            octree.shrink;
            for s = 1:obj.starNumber
                AccVectTemp = [0,0,0];
                for b = 1:octree.binCount
                    switch obj.isFar(octree,b,s)
                        case true
                            binMass = obj.addBinMass(octree,b);
                            binPosDiffVect = obj.binPosVect(octree,b) - obj.positionVectArray(s,:);
                            sjAccVectTemp = obj.G * binMass * binPosDiffVect/(sum(binPosDiffVect.^2))^1.5;
                            AccVectTemp = AccVectTemp + sjAccVectTemp;
                        case false
                            [binStarNumber,binStarMatchIndex] = obj.binCheck(octree,b);
                            for j = 1:binStarNumber
                                if binStarMatchIndex(j) == s
                                    % Do nothing.
                                else
                                    jMass = obj.massStarArray(binStarMatchIndex(j));
                                    sjPosDiffVect = obj.positionVectArray(binStarMatchIndex(j),:) - obj.positionVectArray(s,:);
                                    sjAccVectTemp = obj.G * jMass * sjPosDiffVect/(sum(sjPosDiffVect.^2))^1.5;
                                    AccVectTemp = AccVectTemp + sjAccVectTemp;
                                end
                            end
                    end
                end
                obj.starCluster(s).changeAcceleration(AccVectTemp);
            end
        end
        
        function bool = isFar(obj,octree,b,s)
            if octree.starBinMatch(s) == b
                bool = false;
            else
                binCenter = obj.binPosVect(octree,b);
                sizeComponent = [diff([octree.binBoundary(b,1),octree.binBoundary(b,2)]),...
                                 diff([octree.binBoundary(b,3),octree.binBoundary(b,4)]),...
                                 diff([octree.binBoundary(b,5),octree.binBoundary(b,6)])];
                binSize = max(sizeComponent);
                if binSize^2/sum((binCenter-obj.starCluster(s).position).^2) < obj.farThreshold % Far enough. This is an adjustable threshold value.
                    bool = true;
                else
                    bool = false;
                end
            end
        end
        
        function binCenter = binPosVect(obj,octree,b)
            % USE center of mass
            [binStarNumber,binStarMatchIndex] = obj.binCheck(octree,b);
            binCenterNumerator = [0,0,0];
            for i = 1:binStarNumber
                ss = binStarMatchIndex(i); % ss is the star's index
                binCenterNumerator = binCenterNumerator + obj.massStarArray(ss)*obj.positionVectArray(ss);
            end
            binCenter = binCenterNumerator/obj.addBinMass(octree,b);
        end
        
        function binMass = addBinMass(obj,octree,b)
            [binStarNumber,binStarMatchIndex] = obj.binCheck(octree,b);
            binMass = 0;
            for ss = 1:binStarNumber
                binMass = binMass + obj.massStarArray(binStarMatchIndex(ss));
            end
        end
        
        function [binStarNumber,binStarMatchIndex] = binCheck(~,octree,b)
            binStarMatchCheck = (octree.starBinMatch == b); % vertical array
            binStarNumber = nnz(binStarMatchCheck);
            [binStarMatchIndex,~,~] = find(binStarMatchCheck); % vertical array, binStarNumber elements
        end
    end
end

