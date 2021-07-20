classdef Body < handle
    properties
        mass = [];
        position = [];
        velocity = [];
        acceleration = [];
        accelerationOld = [];
    end
    methods
        function this = Body(initialState)
            if nargin == 0
            else
                this.mass = initialState.mass;
                this.position = initialState.position;
                this.velocity = initialState.velocity;
                this.acceleration = initialState.acceleration;
            end
        end
        function advance(this,dt)
            this.changeVelocity(dt);
            this.changePosition(dt);
        end
        function changePosition(this,dtime)
            %Verlet Algorithm
            this.position = this.position + this.velocity*dtime + 0.5*this.acceleration*dtime^2;
        end
        
        function changeVelocity(this,dtime)
            %Verlet Algorithm
            this.velocity = this.velocity + 0.5*(this.acceleration+this.accelerationOld)*dtime;
        end
    end
end