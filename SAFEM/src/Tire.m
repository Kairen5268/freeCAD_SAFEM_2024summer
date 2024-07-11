classdef Tire
    %TIRE dimension    
    properties
        tire_name
        rib_width
        tire_width
        tire_meshing
        transverse_distribution
        node_width
    end
    
    methods
        function obj = Tire(inputStr,inputVec)
            obj.tire_name = inputStr;
            switch inputStr
                case 'TSD'
                    [obj.rib_width, obj.tire_width, obj.tire_meshing]...
                        = obj.TSD_tire();
                case 'WBT 445'
                    [obj.rib_width, obj.tire_width, obj.tire_meshing,...
                        obj.transverse_distribution,obj.node_width]...
                        = obj.WBT_445();
                otherwise

            end
        end


    end

    methods (Static)

        function [rib_width, tire_width, tire_meshing] = TSD_tire()
            a1 = 264/1000; a2 = 50/1000; a3 = 11/1000; % m, dimension along X
            temp_rib = (a1-4*a3)/5;
            rib_width = temp_rib*ones(1,10); % assume width of each rib is same
            tire_width = a1;
            tire_meshing = [temp_rib/2*ones(1,2) a3 temp_rib/2*ones(1,2) a3...
                            temp_rib/2*ones(1,2) a3 temp_rib/2*ones(1,2) a3...                                                  
                            temp_rib/2*ones(1,2) a2/2*ones(1,2)...                    
                            temp_rib/2*ones(1,2) a3 temp_rib/2*ones(1,2) a3...                     
                            temp_rib/2*ones(1,2) a3 temp_rib/2*ones(1,2) a3...                     
                            temp_rib/2*ones(1,2)];
        end

        % WBT_445 has eight ribs
        function [rib_width, tire_width, tire_meshing, ...
                transverse_distribution,node_width] = WBT_445()
            rib_width = [40 34 40 40 40 40 34 40]/1000;
            rib = rib_width;
            gr=8.4/1000;
            tire_meshing = [rib_width(1)/2*ones(1,2) gr rib_width(2)/2*ones(1,2) gr ...
                rib(3)/2*ones(1,2) gr rib(4)/2*ones(1,2) gr rib(4)/2*ones(1,2) gr ...
                rib(3)/2*ones(1,2) gr rib(2)/2*ones(1,2) gr rib(1)/2*ones(1,2)];
            tire_width = sum(tire_meshing);
            transverse_distribution = [rib_width(1) gr rib_width(2) gr ...
                rib(3) gr rib(4) gr rib(4) gr rib(3) gr rib(2) gr rib(1)];
            node_width = [rib(1)/4*[1 2 1] rib(2)/4*[1 2 1] rib(3)/4*[1 2 1] ...
                rib(4)/4*[1 2 1] rib(4)/4*[1 2 1] rib(3)/4*[1 2 1] rib(2)/4*[1 2 1]...
                rib(1)/4*[1 2 1]];
        end

                
    end    
end

