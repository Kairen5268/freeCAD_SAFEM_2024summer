classdef Dimension
    %% determine the cross-section dimension and layer number
    properties
        layer_number
        layer_thickness
        interface_thickness
        Lx    
    end
    
    methods
        function obj = Dimension(thickness,lx)
            if ~isempty(thickness)
               obj.layer_thickness = thickness;
               obj.interface_thickness = 0.001; % assume as 1 mm
               obj.Lx = lx;
               obj.layer_number=length(thickness);
            end
        end
    end
end