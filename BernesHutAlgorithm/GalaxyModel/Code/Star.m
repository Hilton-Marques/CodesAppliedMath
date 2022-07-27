classdef Star < handle
    % Star consists of properties and methods of a star itself. 
    % Star objects are the fundamental bricks of a star cluster (for now).
    %   Properties: mass, position, velocity, acceleration
    %   Methods:    Star, changeVelocity, changePosition
    %
    %   Verlet Algorithm is applied since V2.0
    %   Verlet Algorithm has good numerical stability and preservation, 
    %     and there is no significant additional computational cost over 
    %     the originally used Euler method. And the global error is of 
    %     order two in stead of one (Wikipedia).
    %   - - - - - - - - - - - -
    %   Author: Haihong
    %           Valentine's Day, 2015
    %   inspired by Salman Khan's videos on astronomy at Khan Academy.
    
    properties
        mass
    end
    properties(Hidden)
        position        % Vector
        velocity        % Vector
        acceleration    % Vector
        accelerationOld % Vector
    end
    
    methods
        function obj = Star(starState) % Constructor
            if nargin == 0
                % Do nothing.
            elseif nargin == 1
                obj.mass = starState.mass;
                obj.position = starState.position;
                obj.velocity = starState.velocity;
                obj.acceleration = starState.acceleration;
            else
                disp('Wrong input number of arguments.\n');
            end
        end
        
        function changePosition(obj,dtime)
            %Verlet Algorithm
            obj.position = obj.position + obj.velocity*dtime + 0.5*obj.acceleration*dtime^2;
        end
        
        function changeVelocity(obj,dtime)
            %Verlet Algorithm
            obj.velocity = obj.velocity + 0.5*(obj.acceleration+obj.accelerationOld)*dtime;
        end
        
        function changeAcceleration(obj,gravityCausedAcceleration)
            obj.accelerationOld = obj.acceleration;
            obj.acceleration = gravityCausedAcceleration;
        end
    end
    
end

