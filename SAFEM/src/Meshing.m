classdef Meshing
    %MESHING matrix Dx and Dy    
    properties
        Dx
        Dy
        layer_element
        element_node
        nodecoordinate
        interface_node
        GDof
        loading_index 
    end
    
    methods
        function obj = Meshing(dimension,tire)
            if isa(dimension, 'Dimension') && isa(tire, 'Tire')
                [obj.Dx, obj.Dy,coordinate_begin] = obj.a_Q4mesh_Seed(dimension.Lx, ...
                    tire.tire_width,...
                    dimension.layer_thickness,tire.tire_meshing);
                [obj.layer_element,obj.element_node,...
                    obj.nodecoordinate,obj.interface_node]=...
                    obj.a_Q4mesh_Meshing(obj.Dx,obj.Dy,...
                    dimension.interface_thickness);
                obj.GDof = 3*length(obj.nodecoordinate);
                obj.loading_index = obj.Loading_index(obj.nodecoordinate,...
                    dimension.layer_thickness,coordinate_begin,...
                    dimension.interface_thickness,tire.tire_meshing,...
                    tire.node_width,tire.tire_name);
            else
                error('Error of Meshing inputs');
            end
            
        end 

        %% loading index and width
        function loading_index = Loading_index(~,nodecoordinate,...
                thickness,coordinate_begin,d_interface,dx_tire,node_width,tire_name)
        
        switch tire_name
            case 'WBT 445'
               for i=1:length(nodecoordinate)
                    if (abs(nodecoordinate(i,1)-coordinate_begin)<=1e-5&&...
                            abs(nodecoordinate(i,2)-sum(thickness)-...
                            d_interface*(length(thickness)-1))<=1e-5)
                        temp1=i;  
                    elseif (abs(nodecoordinate(i,1)-(coordinate_begin+sum(dx_tire)))<=1e-5&&...
                            abs(nodecoordinate(i,2)-sum(thickness)-...
                            d_interface*(length(thickness)-1))<=1e-5) 
                        temp2=i;
                    end
                end
                loading_index=temp1:temp2;
                loading_width=node_width;
                loading_index=[loading_index' loading_width'];
         end
        end

        %% Q4 meshing seed
        function [Dx,Dy,coordinate_begin] = a_Q4mesh_Seed(~,Lx,d_load,thickness,dx_tire)
        N=length(thickness); % layer number
        % default: the load width has 4 element, and the desne field has three load
        % width
        dx_load=[1/4 1/4 1/4 1/8 1/16 1/16]*d_load; 
        Lx_remain=(Lx-d_load*2-sum(dx_tire))/2;
        coordinate_begin=Lx_remain+d_load;
        % default 2 element
        dx1=[1/4 1/4 1/4 1/8 1/8]*Lx_remain;
        % dx1=[1/2*Lx_remain,1/4*Lx_remain,1/4*Lx_remain]; % ,1/8*Lx_remain(1),1/8*Lx_remain(1)];
        % transverse element size
        dx=[dx1,dx_load,dx_tire,fliplr(dx_load),fliplr(dx1)];
        Dx=zeros(N,length(dx));
        for i=1:N
            Dx(i,:)=dx;
        end
        % consider usual conditions 4 5 6 layers
        switch (N)
            case 4
                Dy=zeros(4,10);
                % default 4 element
                Dy(1,1:4)=[1/2*thickness(1),1/4*thickness(1),1/8*thickness(1),1/8*thickness(1)];
                % default 4 element
                Dy(2,1:4)=[1/2*thickness(2),1/4*thickness(2),1/8*thickness(2),1/8*thickness(2)];
                % default 6 element
                Dy(3,1:6)=thickness(3)/6*ones(1,6);
                % default 8 element
                Dy(4,1:10)=[(thickness(4)-0.050)/6*ones(1,6) 0.05/4*ones(1,4)];
            case 5
                Dy=zeros(5,10);
                % default 4 element
                Dy(1,1:4)=[1/2*thickness(1),1/4*thickness(1),1/8*thickness(1),1/8*thickness(1)];
                % default 4 element
                Dy(2,1:4)=[1/4*thickness(2),1/4*thickness(2),1/4*thickness(2),1/4*thickness(2)];
                % default 6 element
                Dy(3,1:6)=thickness(3)/6*ones(1,6);
                % default 8 element
                Dy(4,1:10)=thickness(4)/10*ones(1,10);
                % default 10 element
                Dy(5,1:4)=thickness(5)/4*ones(1,4);
            case 6
                Dy=zeros(6,10);
                % default 4 element
                Dy(1,1:4)=[1/2*thickness(1),1/4*thickness(1),1/8*thickness(1),1/8*thickness(1)];
                % default 4 element
                Dy(2,1:4)=[1/2*thickness(2),1/4*thickness(2),1/8*thickness(2),1/8*thickness(2)];
                % default 4 element
                Dy(3,1:4)=[1/2*thickness(3),1/4*thickness(3),1/8*thickness(3),1/8*thickness(3)];
                % default 6 element
                Dy(4,1:6)=thickness(4)/6*ones(1,6);
                % default 8 element
                Dy(5,1:8)=thickness(5)/8*ones(1,8);
                % default 10 element
                Dy(6,1:10)=thickness(6)/10*ones(1,10);
        end
        end
        %% Q4 meshing
        function [layer_element,element_node,nodecoordinate,interface_node]=...
            a_Q4mesh_Meshing(~,Dx,Dy,d_interface)
        % 2021.3.4 nanjing, SAPAVE 2.0 (finished)
        % generate Q4 meshing & element data from seed function.
        % element size martix Dx(N, ) & Dy(N, ), N is the layer number.
        % meshing generation
        % layer number 
        N=size(Dx,1);
        layer_element=zeros(1,N); % element number of each layer
        Elements=0; Nodes=0;
        interface_node=zeros(N-1,4); % interface_node number
        layer_nodes=zeros(1,N+1);
        for i=1:N
            x=nnz(Dx(i,:)); y=nnz(Dy(i,:)); % the number of nonzero elements
            Elements=Elements+x*y;
            Nodes=Nodes+(x+1)*(y+1);
            layer_element(i)=Elements;
            if i<N
                interface_node(i,:)=[Nodes-x Nodes Nodes+1 Nodes+1+x];
            end
            layer_nodes(i+1)=Nodes;
        end
        % return variable
        element_node=zeros(Elements,4);
        xc=zeros(Nodes,1); yc=zeros(Nodes,1);
        nodecoordinate=zeros(Nodes,2);
        % meshing iteration (direction: down to up)
        Nt=0; Yc=0; 
        layer_element_temp=[0 layer_element];
        for i=1:N
            dx=Dx(i,find(Dx(i,:))); dy=Dy(i,find(Dy(i,:)));  %#ok<*FNDSB>
            nnodex=length(dx)+1; nnodey=length(dy)+1;
            % element-ndoe iteration
            for elerow=1:length(dy)
                for elecol=1:length(dx)
                    % add initial element number for each layer !!!
                    eleno=(elerow-1)*length(dx)+elecol+layer_element_temp(i);
                    % add initial node number for each layer !!!
                    temp=(elerow-1)*length(dx)+elecol+layer_nodes(i);  
                    element_node(eleno,1)=temp+(elerow-1);
                    element_node(eleno,2)=element_node(eleno,1)+1;
                    element_node(eleno,3)=element_node(eleno,2)+nnodex;
                    element_node(eleno,4)=element_node(eleno,1)+nnodex;
                end
            end
            % node coordinate iteration
            if i>1
            temp_interface=d_interface;
            else
                temp_interface=0;
            end
        
            dx_temp=[0 dx]; dy_temp=[0 dy];
            for row=1:nnodey
                for col=1:nnodex 
                    nt=Nt+(row-1)*nnodex+col;
                    xc(nt)=sum(dx_temp(1:col)); yc(nt)=Yc+sum(dy_temp(1:row))+temp_interface;
                end
            end
            Nt=nt; Yc=yc(nt);
        end
        % combining node data into a matrix nodecoordinate=[xc,yc]
        for i=1:Nt
            nodecoordinate(i,1)=xc(i); nodecoordinate(i,2)=yc(i);
        end
        end
    end
end

