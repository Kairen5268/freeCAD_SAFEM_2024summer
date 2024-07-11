classdef Interface
    %INTERFACE Summary of this class goes here  
    
    properties
        stiffness
        strength
    end
    
    methods
        function obj = Interface(meshing,fourier,strength)
            if nargin >= 3 
                obj.strength = strength;  
            else
                temp = length(meshing.layer_element) - 1; % interface num
                obj.strength = [10^7*ones(1,temp-2) 1.0*10^7 1.0e12;% transverse
                                10^7*ones(1,temp-2) 1.0*10^7 1.0e12];% longitudinal
            end           
            obj.stiffness = obj.Interface_spring(...
                 meshing.interface_node,meshing.nodecoordinate,...
                 fourier.Lz);
            obj.stiffness = sparse(obj.stiffness);
        end 
    
        function [interface] = Interface_spring(...
                 obj,interface_node,nodecoordinate,Lz)
            Kx = obj.strength(1,:); Kz = obj.strength(2,:);
            % 2021.3.26 nanjing, SAPAVE 2.0 (finished)
            % generate the martix of spring interface
            % Input: meshing data (interface_node,nodecoordinate,Dx)
            %        interface_node (N-1,4) contains four node number of the interface
            %        Fourier parameters (Lz)
            %        spring parameters (kx,ky,kz) verctor, relating to the number of
            %        interface. kx & kz is the density in the contanct surface !!!
            %        default the Ky as 1e20 to simulate no seperation and no penetration
            % direction: from down to up
            % Note: energy is stored in the dof of martix, for the interface, energy is
            % stored in up-down dof by half.
            % generate the interface martix 
            N=size(interface_node,1); % the number of interface
            Ky=15e12*ones(1,N); % default Ky
            nodes=length(nodecoordinate); Gdof=3*nodes;
            temp=[];
            for i=1:nodes
                if nodecoordinate(i,2)==0
                    temp=[temp nodecoordinate(i,1)]; %#ok<AGROW>
                end
            end
            dx=zeros(1,length(temp)-1);
            for i=1:length(dx)
                dx(i)=temp(i+1)-temp(i);
            end
            dx=[0,dx,0];
            interface=zeros(Gdof);
            % iteration according to interface
            for n=1:N
                a=[interface_node(n,1),interface_node(n,2),interface_node(n,3),interface_node(n,4)];
                node_pair=[a(1):a(2);a(3):a(4)]'; % interface node number
                for i=1:length(node_pair)
                    % 1/2 energy, (dx+dx)/2 linear density
                    tempx=Lz/2*[Kx(n)*(dx(i)+dx(i+1))/2 -Kx(n)*(dx(i)+dx(i+1))/2; 
                                         -Kx(n)*(dx(i)+dx(i+1))/2  Kx(n)*(dx(i)+dx(i+1))/2];
                    interface(node_pair(i,:),node_pair(i,:))=interface(node_pair(i,:),node_pair(i,:))+tempx;      
                    tempy=Lz/2*[Ky(n)*(dx(i)+dx(i+1))/2 -Ky(n)*(dx(i)+dx(i+1))/2; 
                                         -Ky(n)*(dx(i)+dx(i+1))/2   Ky(n)*(dx(i)+dx(i+1))/2];
                    interface(node_pair(i,:)+nodes,node_pair(i,:)+nodes)=...
                        interface(node_pair(i,:)+nodes,node_pair(i,:)+nodes)+tempy;        
                    tempz=Lz/2*[Kz(n)*(dx(i)+dx(i+1))/2 -Kz(n)*(dx(i)+dx(i+1))/2; 
                                         -Kz(n)*(dx(i)+dx(i+1))/2  Kz(n)*(dx(i)+dx(i+1))/2];
                    interface(node_pair(i,:)+2*nodes,node_pair(i,:)+2*nodes)=...
                        interface(node_pair(i,:)+2*nodes,node_pair(i,:)+2*nodes)+tempz; 
                end
            end
        end
    end  
end

