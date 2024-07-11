classdef Temperature
    %TEMPERATURE Summary of this class goes here
    
    properties
        one_dimension
    end
    
    methods
        function obj = Temperature(inputStr,inputArg)
            switch inputStr
                case 'input'
                    obj.one_dimension = Temperature.inter1D(inputArg);




                otherwise

            end
            
        end
    end
        
    methods (Static)
        function distribution = inter1D(temperature_profile)
            depth = 0:0.1:1000; % mm
            Temperature_profile=[0 temperature_profile(1,2); ...
                temperature_profile; 1000 temperature_profile(end,2)]; % mm and ℃
            yi = interp1(Temperature_profile(:,1), ...
                Temperature_profile(:,2), depth,'linear');
            distribution = [depth' yi'];
        end

    end

end

