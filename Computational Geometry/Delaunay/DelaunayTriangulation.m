classdef DelaunayTriangulation < handle
    properties
        triangles
    end
    methods
        function this = DelaunayTriangulation(tri)
            this.triangles = Triangle.empty;
            this.triangles(1) = tri;
        end
        function addPoint(this,p)
            badTriangles = Triangle.empty;
            for j = 1:length(this.triangles)
                T = this.triangles(j);
                if (T.isInCircumcircle(p))
                    badTriangles(end+1) = T;
                end
            end
            %             figure
            %             hold on
            %                 this.plotT(this.triangles,[0,0,1]);
            %                 plot(p(1),p(2),'o','Color','blue');
            %             hold off
            %                         figure
            %                         hold on
            %                         this.plotT(this.triangles,[0,0,1]);
            %                         this.plotT(badTriangles,[1,0,0],flag);
            %                         hold off
            
            boundary = this.getBoundary(badTriangles);
            %                         figure
            %                         hold on
            %                         this.plotT(this.triangles,[0,0,1]);
            %                         this.plotT(badTriangles,[1,0,0],flag);
            %                         this.plotE(boundary);
            %                         hold off
            indexBad = ~ismember(this.triangles , badTriangles);
            this.triangles = this.triangles(indexBad);
            newTriangles = Triangle.empty;
            for i = 1:length(boundary)
                edge = boundary(i);
                newTri = Triangle(p,edge.p1,edge.p2);
                newTri.neigh(1) = edge.T;
                newTri.edges(1).T = edge.T;
                newTriangles(end+1) = newTri;
            end
            N = length(newTriangles);
            for i = 1:N
                newTriangles(i).neigh(3) = ...
                    newTriangles(mod(mod(i-2,N) + N,N)+1);
                newTriangles(i).edges(3).T = newTriangles(i).neigh(3);
                newTriangles(i).neigh(2) = ...
                    newTriangles(mod(i,N)+1);
                newTriangles(i).edges(2).T = newTriangles(i).neigh(2);
            end
            
            
            % Update edges of Old triangles. This could be avoid with a
            % shared pointer
            this.updateEdges(newTriangles);
            this.triangles(end+1:end+N) =  newTriangles;
            figure
            hold on
            this.plotT(this.triangles,[0,0,1]);
            plot(p(1),p(2),'o','Color','blue');
            hold off
        end
        function boundary = getBoundary(~,badTriangles)
            boundary = Edge.empty;
            if isempty(badTriangles)
                return;
            end
            T = badTriangles(1);
            edge = 1;
            while (true)
                if (length(boundary) > 1)
                    if (boundary(1) == boundary(end))
                        break
                    end
                end
                index = ismember(badTriangles,T.neigh(edge));
                if (~index)
                    boundary(end+1) = T.edges(edge);
                    edge = mod((edge),3) + 1;
                else
                    last = T;
                    T = badTriangles(index);
                    [~,pos] = ismember(last,T.neigh);
                    edge = mod((pos),3) + 1;
                end
            end
            boundary(end) = [];
        end
        function updateEdges(this,newTri)
            allEdges = [newTri(:).edges];
            if isempty(allEdges)
                return
            end
            edgesP1 = [allEdges.p1]';
            edgesP2 = [allEdges.p2]';
            allEdgesp1 = [edgesP1(1:2:end),edgesP1(2:2:end)];
            allEdgesp2 = [edgesP2(1:2:end),edgesP2(2:2:end)];
            allEdgesP = [allEdgesp1,allEdgesp2];
            for i = 1:length(this.triangles)
                T = this.triangles(i);
                for j = 1:3
                    edge = T.edges(j);
                    [~,pos1] = ismember(allEdgesP,[edge.p1,edge.p2],'rows');
                    [~,pos2] = ismember([edge.p2,edge.p1],allEdgesP,'rows');
                    if pos1 > 0
                        T.neigh(j) = newTri(fix((pos1-1)/3) + 1);
                        T.edges(j).T = T.neigh(j);
                        continue;
                    elseif pos2 > 0
                        T.neigh(j) = newTri(fix((pos2-1)/3) + 1);
                        T.edges(j).T = T.neigh(j);
                    end
                end
            end
        end
        function plotT(~,triangles,c,flag)
            for i = 1: length(triangles)
                triangles(i).plot(c);
                if nargin > 3
                    pMed = (triangles(i).points(1,:) + ...
                        triangles(i).points(2,:) + ...
                        triangles(i).points(3,:))/3;
                    %text(pMed(1),pMed(2),num2str(i));
                end
            end
        end
        function plotE(~,boundary)
            for i = 1: length(boundary)
                boundary(i).plot;
                pMed = (boundary(i).p1 + boundary(i).p2)/2;
                text(pMed(1),pMed(2),num2str(i));
            end
        end
        function plotNeigh(~,bd)
            for i = 1:3
                neig = bd.neigh(i);
                if ~isempty(neig.edges)
                    neig.plot([1,0,1]);
                end
            end
        end
    end
end
