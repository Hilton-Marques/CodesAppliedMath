classdef Point < handle
    properties
        coord
        id
        elements = Face.empty;
        idTemp;
        u = [];
        bd = false;
        markTemp = false;
        closeElements = Face.empty;
        star = Point.empty;
        hedStar = Hed.empty;
        hed0
        new_coord = [];
    end
    methods
        function this = Point(coord,id)
            if nargin == 2
                this.coord = coord;
                this.id = id;
            end
        end
        function star = getStar(this,solid)
            star = [];
            v0 = this.hed0.inc(2); 
            v1 = -1;
            hedi = this.hed0; 
            star(end+1) = v0;
            while (v1 ~= v0)                
                twin = solid.edges(hedi.edgeId).getTwin(hedi.id);
                elementi = (twin.elId);
                next = solid.heds(twin.heNext);
                v1 = next.inc(2);
                hedi = next;
                star(end+1) = v1;
                %elements(end+1) = elementi;
            end
            star(end) = [];
            %this.elements = elements;
        end
        function out = getSolidAngle(this)
            out = 0;
            size = length(this.elements);
            for element = this.elements 
                element.markTemp = true;
            end
            for elementi = this.elements                
                ni = elementi.geometry.getNormal(0,0);
                for elementj = elementi.adjacentElementsEdges
                    if elementj.markTemp
                        nj = elementj.geometry.getNormal(0,0);
                        angle = this.getAngle(ni,-nj);
                        out = out + angle;
                    end
                end
                elementi.markTemp = false;
            end
             out = (out - pi)/(4*pi);
            if (size ~= 3)
              out = out - 0.75;
            end
        end
        
        function [im,rot] = turnAround(this,solid, rot , az)
            heds = solid.heds;
            points = solid.points;
            elements = solid.elements;
            edges = solid.edges;
            this.plot();
            count = 0;
            s_count = int2str(count);   
            rot = rot + 3;
            view(rot,az);
            drawnow;
            countFig = 1;
            frame = getframe(gca);
            im{countFig} = frame2im(frame);
            countFig = countFig + 1;
            
            %exportgraphics(gca,strcat('Neigh-point',s_count ,'.png'),'Resolution',1000);
            
            %find fist hed
            flag = false;
            for (el = this.elements)
                for (hed = el.heds)
                    if (hed.inc(1) == this.id)
                        hed0 = hed;
                        flag = true;
                        break;
                    end
                end
                if (flag)
                    break;
                end
            end
            
            v0 = hed0.inc(2); 
            v1 = -1;
            %hed0.plot(heds,points,'red');
            frame = getframe(gca);
            im{countFig} = frame2im(frame);
            countFig = countFig + 1;
            count = count + 1;
            s_count = int2str(count);
            %exportgraphics(gca,strcat('Neigh-point',s_count ,'.png'),'Resolution',1000);
            elementi = elements(hed0.elId);          
            elementi.plot('b',0.6,true);
            rot = rot + 3;
            view(rot,az);
            drawnow;
            frame = getframe(gca);
            im{countFig} = frame2im(frame);
            countFig = countFig + 1;
            hedi = hed0; 
            while (v1 ~= v0)
                %elementi = elements(hedi.elId);
                %elementi.plot('r',0.8,true);
                twin = edges(hedi.edgeId).getTwin(hedi.id);
                %twin.plot(heds,points,'green');
                rot = rot + 3;
                view(rot,az);
                drawnow;
                frame = getframe(gca);
                im{countFig} = frame2im(frame);
                countFig = countFig + 1;
                elementi = elements(twin.elId);
                elementi.plot('b',0.6,true);
                frame = getframe(gca);
                im{countFig} = frame2im(frame);
                
                countFig = countFig + 1;
                count = count + 1;
                s_count = int2str(count);
                %exportgraphics(gca,strcat('Neigh-point',s_count ,'.png'),'Resolution',1000);
                next = heds(twin.heNext);
                %next.plot(heds,points,'red');
                rot = rot + 3;
                view(rot,az);
                drawnow;                
                frame = getframe(gca);
                im{countFig} = frame2im(frame);
                countFig = countFig + 1;
                count = count + 1;
                s_count = int2str(count);
                %exportgraphics(gca,strcat('Neigh-point',s_count ,'.png'),'Resolution',1000);
                v1 = next.inc(2);
                hedi = next;
            end
            
        end
        function h = plot(this)
            p0 = this.coord;
            h = plot3(p0(1),p0(2),p0(3),'o','MarkerSize',5,'Color','cyan','MarkerFaceColor','cyan');
        end
    end
        methods (Static)
        function out = getAngle(v1,v2)
            v1 = v1 / norm(v1);
            v2 = v2/ norm(v2) ;
            out = acos(dot(v1,v2));
        end
    end
end