classdef ResultsFigure
    %RESULTSFIGURE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = ResultsFigure(inputStr,results,meshing)
            if contains(inputStr, 'cross-section', 'IgnoreCase', true)
               ResultsFigure.cross_sectionfigure(inputStr,results,meshing);
            end
        end
        
    end

    methods(Static)

        function damagefigure(mapping, damage)
            fig = figure("Name",'Damage Distribution');
            fig.Position = [100, 100, 600, 400];
            XX=[]; YY=[];
            for i = 1:length(mapping.element_node)
                indicate = mapping.element_node(i,:);
                XX0 = mapping.nodecoordinate(indicate,1);
                YY0 = mapping.nodecoordinate(indicate,2);
                XX = [XX XX0]; YY = [YY YY0];
            end
            XX=XX-0.5*max(XX(:)); YY=YY-max(YY(:));
            patch(XX,YY,damage');
            colorbar;
            caxis([0 1]);
            xlabel('Coordinate X (m)','Fontweight','bold'); 
            ylabel('Coordinate Y (m)','Fontweight','bold'); 
            title('Damage Distribution')
            colormap('turbo');
        end

        function cross_sectionfigure(inputStr,results,meshing)
            fig = figure("Name",inputStr);
            fig.Position = [100, 100, ...
                800*[1, max(meshing.nodecoordinate(:,2))...
                /max(meshing.nodecoordinate(:,1))]];
            X=[]; Y=[];
            for i = 1:length(meshing.element_node)
                indicate = meshing.element_node(i,:);
                X0 = meshing.nodecoordinate(indicate,1);
                Y0 = meshing.nodecoordinate(indicate,2);
                X = [X X0]; Y = [Y Y0];
            end
            X=X-0.5*max(X(:)); Y=Y-max(Y(:));
            if contains(inputStr, 'transverse shear', 'IgnoreCase', true)
                index = 4;
            elseif contains(inputStr, 'longitudinal shear', 'IgnoreCase', true)
                index = 5;
            elseif contains(inputStr, 'longitudinal', 'IgnoreCase', true)
                index = 3;
            elseif contains(inputStr, 'transverse', 'IgnoreCase', true)
                index = 1;
            elseif contains(inputStr, 'compressive', 'IgnoreCase', true)
                index = 2;
            end
            if contains(inputStr, 'strain', 'IgnoreCase', true)
                response_temp = results.strain_fixed(:,:,index,:)*1e6;
            elseif contains(inputStr, 'stress', 'IgnoreCase', true)
                response_temp = results.stress_fixed(:,:,index,:);
            end
            [maxstress,position] = max(abs(response_temp(:)));
            [dim1,dim2,dim3,dim4] = ind2sub(size(response_temp),position);  
            patch(X,Y,squeeze(response_temp(:,:,1,dim4))');
            colorbar;
            xlabel('Coordinate X (m)','Fontweight','bold'); 
            ylabel('Coordinate Y (m)','Fontweight','bold'); 
            title([inputStr, ' (micro)'])
            colormap('turbo'); 
            hold on
            scatter(meshing.nodecoordinate(meshing.element_node(dim1,dim2),1)-...
                0.5*max(meshing.nodecoordinate(:,1)),...
                meshing.nodecoordinate(meshing.element_node(dim1,dim2),2)-...
                max(meshing.nodecoordinate(:,2)),'blue','LineWidth',1);
        end

        function Mappingcomparision(mapping,meshing,ylimit,xlimit)
            if nargin == 3
                yl = ylimit;
            end
            if nargin == 4
                yl = ylimit;
                xl = xlimit;
            end
            fig = figure("Name",'Mappingcomparision');
            fig.Position = [100, 100, 1200, 400];
            subplot(1,2,1)
            X=[]; Y=[];
            for i = 1:length(meshing.element_node)
                indicate = meshing.element_node(i,:);
                X0 = meshing.nodecoordinate(indicate,1);
                Y0 = meshing.nodecoordinate(indicate,2);
                X = [X X0]; Y = [Y Y0];
            end
            X=X-0.5*max(X(:)); Y=Y-max(Y(:));
            patch(X,Y,mapping.selected_results');
            colorbar;
            xlabel('Coordinate X (m)','Fontweight','bold'); 
            ylabel('Coordinate Y (m)','Fontweight','bold'); 
            title('Original Meshing')
            colormap('turbo');
            xlim([-xl xl]);
            ylim([yl 0]);

            subplot(1,2,2)
            XX=[]; YY=[];
            for i = 1:length(mapping.element_node)
                indicate = mapping.element_node(i,:);
                XX0 = mapping.nodecoordinate(indicate,1);
                YY0 = mapping.nodecoordinate(indicate,2);
                XX = [XX XX0]; YY = [YY YY0];
            end
            XX=XX-0.5*max(XX(:)); YY=YY-max(YY(:));
            patch(XX,YY,mapping.transform_results');
            colorbar;
            xlabel('Coordinate X (m)','Fontweight','bold'); 
            ylabel('Coordinate Y (m)','Fontweight','bold'); 
            title('Projection Meshing')
            colormap('turbo');
            xlim([-xl xl]);
            ylim([yl 0]);
        end

        function directpatch(X, Y, Z)
            patch(X,Y,Z);
        end
    end
end

