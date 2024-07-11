classdef Process
    %PROCESS Summary of this class goes here
    
    properties
        speed
        step
        time
        loading_factor
        z_coordinate
    end
    
    methods
        function obj = Process(inputArg1,inputArg2,inputArg3,initial_z,inputArg4)
            %PROCESS Construct an instance of this class
            obj.speed = [zeros(1,10) inputArg1*ones(1,inputArg2+10)];
            obj.step = inputArg2+20;
            obj.time = inputArg3;
            if nargin >= 5 
                obj.loading_factor = inputArg4;  
            else
                obj.loading_factor = [(0.1:0.1:1) ones(1,obj.step-10)]; 
            end
            obj.z_coordinate = zeros(1, obj.step);
            for i =1:obj.step
                obj.z_coordinate(i) = initial_z + sum(obj.speed(1:i)*obj.time);
            end
        end
        
        
    end
end

