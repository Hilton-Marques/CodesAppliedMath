classdef LogClass < handle
    % This is a log class
    %   Properties: logname, fid, starNumber
    %   Private Methods: LogClass
    %   Public Methods:  getInstance, deleteObj, printTime, printContent
    %   - - - - - - - - - - - -
    %   Author: Haihong
    %           Valentine's Day, 2015
    %   inspired by Salman Khan's videos on astronomy at Khan Academy.
    
    properties
        logName    % A string
    end
    properties(Hidden)
        fid        % A file identifier
        starNumber % Used to create a new line
    end
    
    methods(Access = private)
        function obj = LogClass(starNumber) % Constructor, private method
            datestrout = datestr(now,30);
            currentDir = pwd();
            obj.logName = strcat(currentDir,'\Log date',datestrout(1:8),' time',datestrout(10:15),'.txt');
            obj.fid = fopen(obj.logName,'a'); % Open or create for reading and writing, append data to end of file
            obj.starNumber = starNumber;
        end
    end
    
    methods(Static)
        function obj = getInstance(starNumber) % Accepts external call
            persistent localObj;
            if isempty(localObj)|| ~isvalid(localObj)
                localObj = LogClass(starNumber); % If obj does not exist, then create one.
            end
            obj = localObj;
        end
    end
    
    methods
        function deleteObj(obj) % Encapsulate FCLOSE
            fclose(obj.fid);
        end
        
        function printTime(obj,t) % Encapsulate FPRINTF
            fprintf(obj.fid,'%.2f%%',t);
        end
        
        function printContent(obj,s,star) % Encapsulate FPRINTF
            fprintf(obj.fid,'%d\tr\t%f\t%f\t%f\tv\t%f\t%f\t%f\t',...
                               s,...
                               star.position(1),star.position(2),star.position(3),...
                               star.velocity(1),star.velocity(2),star.velocity(3));
            if s == obj.starNumber
               fprintf(obj.fid,'%s\r\n',datestr(now,0)); % Microsoft Notepad requires '\r\n' instead of '\n'.
            end
        end
    end
    
   
end

