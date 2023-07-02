classdef Leaf < handle
    properties
        s;
        nodeCE = [];
        bodies;
        capacity = 0;
        totalMass = 0;
        leafMassXCM = zeros(1,2);
        isDivided = false;
        nBodies = 0;
        children = cell(4,1);
    end
    methods
        function this = Leaf(s,cellCentroid,capacity)
            this.s = s;
            this.nodeCE = cellCentroid;
            this.capacity = capacity; 
            %PreAllocate bodies
            bodiesCreation(capacity) = Body();
            this.bodies = bodiesCreation;
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
            this.children{1,1} = Leaf(s/2,oldCE + [-s/4,-s/4],capacity);
            this.children{2,1} = Leaf(s/2,oldCE +  [s/4,-s/4],capacity);
            this.children{3,1} = Leaf(s/2,oldCE +  [-s/4,s/4],capacity);
            this.children{4,1} = Leaf(s/2,oldCE +  [s/4,s/4] ,capacity);
            if(nBody1 > 0)
                this.children{1,1}.bodies(1:nBody1,1) = [this.bodies(idx1quadrant)];
                this.children{1,1}.updateTotal(nBody1);
            end
            if (nBody2 > 0)
                this.children{2,1}.bodies(1:nBody2,1) = [this.bodies(idx2quadrant)];
                this.children{2,1}.updateTotal(nBody2);
            end
            if (nBody3 > 0)
                this.children{3,1}.bodies(1:nBody3,1) = [this.bodies(idx3quadrant)];
                this.children{3,1}.updateTotal(nBody3);
            end
            if (nBody4 > 0)
                this.children{4,1}.bodies(1:nBody4,1) = [this.bodies(idx4quadrant)];
                this.children{4,1}.updateTotal(nBody4);
            end
        end
        function output = getQuadrant(this,body)
            position = body.position;
            gtMidMask = position >= this.nodeCE;
            output = gtMidMask(:,1) + 2*gtMidMask(:,2) + 1;
        end
        function update(this,body)
            this.nBodies = this.nBodies + 1;
            this.bodies(this.nBodies,1) = body;
            this.leafMassXCM = (this.leafMassXCM*this.totalMass + ...
                body.mass * body.position)/(this.totalMass + body.mass);
            this.totalMass = this.totalMass + body.mass;   
        end
        function output = distSqr(this,body)
            leafCM = this.leafMassXCM;
            r = (leafCM - body.position);
            output = dot(r,r);
        end
        function updateAccelerationOn(this,body,rSqr)
            %Case of the body is much closer or equal to Nodes's centroid
            if (rSqr < 0.002)
                return;
            end
            M = this.totalMass;
            mXP = this.leafMassXCM;
            invDCube = rSqr^(-1.5);
            incrementalAcce = M*(mXP - body.position) * invDCube;
            body.acceleration = body.acceleration + incrementalAcce;
        end
        function updateTotal(this,nBodies)
            this.nBodies = nBodies;
            positions = [this.bodies(:).position];
            bodiesMass = [this.bodies(:).mass]';
            positionsVec = [ [positions(1:2:end)]' [positions(2:2:end)]'];
            this.totalMass = sum(bodiesMass);
            this.leafMassXCM = sum(bodiesMass .* positionsVec ,1) / this.totalMass ;
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
        function [h1,h2] = showCurrent(this,body)
            left_corner_x = this.nodeCE(1,1) - this.s*0.5;
            left_corner_y = this.nodeCE(1,2) - this.s*0.5;
            h1 = rectangle('Position',[left_corner_x,left_corner_y,this.s,this.s],...
                'EdgeColor','g');
            h2 = plot(this.leafMassXCM(1,1),this.leafMassXCM(1,2),'o','Color','red');
            plot(body.position(1,1),body.position(1,2),'o','Color','y');
        end
    end
end