classdef Boundary
    %BOUNDARY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Bottom
        Vertical
    end
    
    methods
        function obj = Boundary(material,meshing,fourier)
            % bottom condition (spring modulus linear density)
            obj.Bottom = Boundary.bottomcondition(...
                material.E(1),material.v(1),meshing.Dx(1,:),...
                3*length(meshing.nodecoordinate),length(meshing.nodecoordinate),...
                meshing.nodecoordinate,fourier.Lz);
            obj.Bottom = sparse(obj.Bottom);
            % vertical condition (damping factor linear density)
            % derived from different layers
            obj.Vertical = Boundary.verticalcondition([material.E sum(material.Ei,2)'],...
                material.v,material.rho,...
                3*length(meshing.nodecoordinate),meshing.nodecoordinate,...
                fourier.Lz,sum(meshing.Dx(1,:)),meshing.Dy);
            obj.Vertical = sparse(obj.Vertical); 
        end
    end
        
    methods (Static)
        function [bottom] = bottomcondition(E,v,Dx,GDof,...
                numbernodes,nodecoordinates,Lz)
            Kv=2.35*E*(Dx/Lz).^0.5/(1-v(1)^2)/2;
            Kh=4.6*E*(Dx/Lz).^0.5/((1+v(1))*(2-v(1)))/2;
            bottom=zeros(GDof);
            indice=find(nodecoordinates(:,2)==0);
            Kv_left=[0 Kv]; Kv_right=[Kv 0]; 
            temp_Kv=(Kv_left+Kv_right)*Lz/2;
            Kh_left=[0 Kh]; Kh_right=[Kh 0]; 
            temp_Kh=(Kh_left+Kh_right)*Lz/2;
            for i=1:length(indice)
                bottom(indice(i)+numbernodes,indice(i)+numbernodes)=...
                    temp_Kv(i);
                bottom(indice(i),indice(i))=temp_Kh(i);
                bottom(indice(i)+2*numbernodes,indice(i)+2*numbernodes)=...
                    temp_Kh(i);
            end
        end

        function [vertical] = verticalcondition(E,v,rho,GDof,...
                nodecoordinates,Lz,Lx,Dy)
            numbernodes=length(nodecoordinates);
            vertical=zeros(GDof);
            indice1=find(nodecoordinates(:,1)==0);
            indice2=find(nodecoordinates(:,1)==Lx);
            segment=sum(Dy~=0,2);
            N=length(segment);
            node=1;
            for i=1:N
                Vp=(E(i)*(1-v(i))/((1+v(i))*(1-2*v(i))*rho(i)))^0.5;
                Vs=(E(i)/(2*(1+v(i))*rho(i)))^0.5;
                vertical(indice1(node:node+segment(i)),indice1(node:node+segment(i)))=...
                    diag(Lz*([0 Dy(i,1:segment(i))]+[Dy(i,1:segment(i)) 0])/2*rho(i)*Vp/2);
                vertical(indice1(node:node+segment(i))+numbernodes,indice1(node:node+segment(i))+numbernodes)=...
                    diag(Lz*([0 Dy(i,1:segment(i))]+[Dy(i,1:segment(i)) 0])/2*rho(i)*Vs/2);
                vertical(indice1(node:node+segment(i))+2*numbernodes,indice1(node:node+segment(i))+2*numbernodes)=...
                    diag(Lz*([0 Dy(i,1:segment(i))]+[Dy(i,1:segment(i)) 0])/2*rho(i)*Vs/2);
                vertical(indice2(node:node+segment(i)),indice2(node:node+segment(i)))=...
                    diag(Lz*([0 Dy(i,1:segment(i))]+[Dy(i,1:segment(i)) 0])/2*rho(i)*Vp/2);
                vertical(indice2(node:node+segment(i))+numbernodes,indice2(node:node+segment(i))+numbernodes)=...
                    diag(Lz*([0 Dy(i,1:segment(i))]+[Dy(i,1:segment(i)) 0])/2*rho(i)*Vs/2);
                vertical(indice2(node:node+segment(i))+2*numbernodes,indice2(node:node+segment(i))+2*numbernodes)=...
                    diag(Lz*([0 Dy(i,1:segment(i))]+[Dy(i,1:segment(i)) 0])/2*rho(i)*Vs/2);
                node=node+segment(i)+1;
            end
        end
    end
end

