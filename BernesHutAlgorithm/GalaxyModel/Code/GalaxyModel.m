classdef GalaxyModel < handle
    % This is the class for the simulator.
    %   Properties:
    %       starNumber, initialState, starCluster, log
    %   Public method (Usage):
    %       create:   g = GalaxyModel(initialState);
    %       simulate: g.simulate(time,dtime);
    %   Protected method:
    %       gravityAffect: classic Newtonian Mechanics is applied;
    %       viewStarCluster: plots starCluster, returning a snapshot
    %
    %   IMPROVEMENT HISTORY:
    %   Verlet Algorithm is applied since Version 2.0
    %   Barnes-Hut Algorithem (Octree Method) is applied since Version 3.0
    %   Gravitational calculation is encapsulated since Version 3.0
    %   Visualization is encapsulated since Version 4.0
    %
    %   PROBLEM OF SINGULARITY:
    %   NOTE that I have not solved the singularity problem of collision.
    %     In this current vision, stars won't collide: they just 'pass
    %     through' each other (with extrememly high speed, because they are
    %     extremely close to each other prior to collision), which is
    %     against physical reality. I originally thought the solution was
    %     easy, but later I found it very tricky.
    %     From literature I have reviewed, the solutions are complicated
    %     mathematically and computationally.
    %   - - - - - - - - - - - -
    %   Author: Haihong
    %           Valentine's Day, 2015
    %   inspired by Salman Khan's videos on astronomy at Khan Academy.
    %   Acknowledgement:
    %   Sven's work Octree is modified and integrated in this program.
    
    properties
        starNumber
        initialState % Array of structs,
        % struct members: mass, position, velocity, acceleration
        starCluster  % Array of Star objects, composition relationship
        log          % A LogClass object, composition relationship
        viewer       % A Viewer object, composition relationship
    end
    
    methods
        function obj = GalaxyModel(initialState) % Constructor
            obj.starNumber = length(initialState);
            obj.initialState = initialState;
            starClusterCreate(obj.starNumber) = Star(); % Just preallocate
            obj.starCluster = starClusterCreate;
            obj.log = LogClass.getInstance(obj.starNumber); % Create a log file
            for i = 1:obj.starNumber
                % Note that starCluster(i) is an object
                obj.starCluster(i) = Star(initialState(i));
            end
            
            % Record to log
            for s = 1:obj.starNumber
                obj.log.printTime(0);
                obj.log.printContent(s,obj.starCluster(s));
            end
            % Constrct a viewer
            obj.viewer = Viewer(obj.starNumber,obj.starCluster,0);
        end
        
        function simulate(obj,time,dtime)
            % Simulate the evolution of the celestial system
            % time: total time
            % dtime: step length
            if ~isa(obj,'GalaxyModel')
                dsp('The input object is not a GalaxyModel object.\n');
                return
            elseif ~isa(obj.log,'LogClass')
                dsp('The property log is not a LogClass object.\n');
                return
            elseif time/dtime ~= fix(time/dtime)
                dsp('Total steps defined by given time and dtime is not an integer.\n');
                return
            elseif obj.starNumber == 0
                dsp('Total number of stars is zero.\n');
                return
            end
            
            % Simualte over time, and record data to log.
            % Note that log is a LogClass object.
            for t = dtime:dtime:time
                obj.gravityAffect; % Accelerations changed by gravity
                for s = 1:obj.starNumber
                    obj.starCluster(s).changeVelocity(dtime); % Second-order Verlet Algorithm applied
                    obj.starCluster(s).changePosition(dtime); % Second-order Verlet Algorithm applied
                    % Record to log at fixed times
                    if (time/dtime <= 200) || ...
                            (time/dtime > 200 && t/(5*dtime) == fix(t/(5*dtime)))
                        % Not print every loop of t, to avoid redundancy.
                        obj.log.printTime(t);
                        obj.log.printContent(s,obj.starCluster(s));
                    end
                end
                obj.viewStarCluster(t);
            end % Simulation done
            obj.log.deleteObj(); % Close the log file
        end
    end
    
    methods(Access = protected)
        function gravityAffect(obj)
            % The law of gravity is used in this function, returning each
            %   star's acceleration (vector). Here only Newtonian Machanics
            %   is applied. No General Relativity, though I have considered
            %   using it.
            % Note that bsxfun(FUNC,A,B) applies the element-by-element
            %   binary operation specified by the function handle FUNC to
            %   matrices A and B. It is very useful in vectorization.
            CONST_G = 1; % Constant of gravitation
            ce = GravityCalculationEngine(CONST_G,obj.starCluster,obj.starNumber);
            obj.starCluster = ce.calculate;
        end
        
        function viewStarCluster(obj,t)
            obj.viewer.update(obj.starCluster,t);
            pause(0.001); % pause for 0.001 s, just for debugging
        end
    end
    
end

