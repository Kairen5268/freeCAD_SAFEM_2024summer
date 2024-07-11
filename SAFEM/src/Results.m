classdef Results
    %RESULTS Summary of this class goes here 
    
    properties
        results_type
        displacement_fixed % fixed location
        velocity_fixed
        acceleration_fixed
        strain_fixed
        stress_fixed
    end
    
    methods
        function obj = Results(inputStr,coordinatez,meshing,process,...
                fourier,solution,material)
            %RESULTS Construct an instance of this class
            obj.results_type = inputStr;
            switch inputStr
                case 'fixed_location'
                    [obj.displacement_fixed,obj.velocity_fixed,...
                        obj.acceleration_fixed] = obj.output_response(...
                        meshing.GDof,process.step,coordinatez,...
                        fourier.Lz,fourier.L,...
                        solution.Displacement,solution.Velocity,...
                        solution.Acceleration);                    
                    straintemp = zeros(size(meshing.element_node,1),...
                        size(meshing.element_node,2),6,process.step);
                        parfor i = 1:process.step
                            straintemp(:,:,:,i) = obj.strain3D(...
                                meshing.element_node,meshing.nodecoordinate,...
                                solution.Displacement(:,i,:),fourier.L,...
                                coordinatez,fourier.Lz);   
                        end
                    obj.strain_fixed = straintemp;
                    stress_elastic = obj.stress3D_elastic(obj.strain_fixed,...
                            meshing.element_node,material.D(:,:,length(material.E)),...
                            [0 meshing.layer_element],process.step);
                    stress_visco = obj.viscoelasticstress3D(obj.strain_fixed,...
                                [0 meshing.layer_element],material.Ei,material.T,...
                                material.D(:,:,length(material.E)+1:end),...
                                material.v,process.time,process.step,...
                                meshing.element_node);
                    obj.stress_fixed = stress_elastic + stress_visco;
            end
        end
        
        function [displacement,velocity,acceleration] = output_response(...
                ~,GDof,n,z,Lz,L,Displacement,Velocity,Acceleration)
            % output displacement, velocity, and acceleration
            displacement=zeros(GDof,n);
            velocity=zeros(GDof,n);
            acceleration=zeros(GDof,n);
            numbernodes=GDof/3;
            for i=1:length(L)
                l=L(i);
            displacement(1:numbernodes,:)=displacement(1:numbernodes,:)+...
                Displacement(1:numbernodes,:,i)*sin(l*pi*z/Lz);
            displacement(1+numbernodes:2*numbernodes,:)=...
                displacement(1+numbernodes:2*numbernodes,:)+...
                Displacement(1+numbernodes:2*numbernodes,:,i)*sin(l*pi*z/Lz); 
            displacement(1+2*numbernodes:3*numbernodes,:)=...
                displacement(1+2*numbernodes:3*numbernodes,:)+...
                Displacement(1+2*numbernodes:3*numbernodes,:,i)*cos(l*pi*z/Lz);
            velocity(1:numbernodes,:)=velocity(1:numbernodes,:)+...
                Velocity(1:numbernodes,:,i)*sin(l*pi*z/Lz);
            velocity(1+numbernodes:2*numbernodes,:)=...
                velocity(1+numbernodes:2*numbernodes,:)+...
                Velocity(1+numbernodes:2*numbernodes,:,i)*sin(l*pi*z/Lz); 
            velocity(1+2*numbernodes:3*numbernodes,:)=...
                velocity(1+2*numbernodes:3*numbernodes,:)+...
                Velocity(1+2*numbernodes:3*numbernodes,:,i)*cos(l*pi*z/Lz);
            acceleration(1:numbernodes,:)=acceleration(1:numbernodes,:)+...
                Acceleration(1:numbernodes,:,i)*sin(l*pi*z/Lz);
            acceleration(1+numbernodes:2*numbernodes,:)=...
                acceleration(1+numbernodes:2*numbernodes,:)+...
                Acceleration(1+numbernodes:2*numbernodes,:,i)*sin(l*pi*z/Lz); 
            acceleration(1+2*numbernodes:3*numbernodes,:)=...
                acceleration(1+2*numbernodes:3*numbernodes,:)+...
                Acceleration(1+2*numbernodes:3*numbernodes,:,i)*cos(l*pi*z/Lz);
            end
        end

        function [ strain ] = strain3D(~,elementnodes,...
                nodecoordinates,displacements,L,z,Lz)
            numbernodes = max(elementnodes(:));
            element = 1:length(elementnodes);
            [gaussweights,gausslocations] = Results.GaussQuadrature();
            numberelements = length(element);
            strain = zeros(numberelements,size(elementnodes,2),6);
            strain0 = zeros(numberelements,size(elementnodes,2),6);
            for j=1:numberelements
                e=element(j);
                indice=elementnodes(e,:);
                elementdof=[indice indice+numbernodes indice+2*numbernodes];
                nn=length(indice);
                for q=1:size(gaussweights,1)
                        pt=gausslocations(q,:);
                        wt=gaussweights(q);
                        xi=pt(1);
                        eta=pt(2);
                    [shapefunction,naturalderivatives] = ...
                        Results.ShapeFunctionQ4(xi,eta);
                    [~,XYderivatives] = Results.Jacobian(...
                        nodecoordinates(indice,:),naturalderivatives);
                    for i=1:length(L)
                        l=L(i);
                        gammal=l*pi/Lz;
                        B=zeros(6,3*nn);
                        B(1,1:nn)=XYderivatives(:,1)'*sin(l*pi*z/Lz);
                        B(2,nn+1:2*nn)=XYderivatives(:,2)'*sin(l*pi*z/Lz);
                        B(3,2*nn+1:3*nn)=-shapefunction*gammal*sin(l*pi*z/Lz);
                        B(4,1:nn)=XYderivatives(:,2)'*sin(l*pi*z/Lz);
                        B(4,nn+1:2*nn)=XYderivatives(:,1)'*sin(l*pi*z/Lz);
                        B(5,nn+1:2*nn)=shapefunction*gammal*cos(l*pi*z/Lz);
                        B(5,2*nn+1:3*nn)=XYderivatives(:,2)'*cos(l*pi*z/Lz);
                        B(6,1:nn)=shapefunction*gammal*cos(l*pi*z/Lz);
                        B(6,2*nn+1:3*nn)=XYderivatives(:,1)'*cos(l*pi*z/Lz);
                        strain0(j,q,:)=B*displacements(elementdof,i);
                        strain(j,q,:)=strain(j,q,:)+strain0(j,q,:);
                    end
                end
            end
        end

        function [ stress ] = stress3D_elastic(~,strain,elementnodes,...
                D,layernumber,n)
            stress = zeros(size(elementnodes,1),size(elementnodes,2),6,n);
            for i = 1:n % analysis step interation
                for j = 1:size(D,3) % D interation
                    Dj = squeeze(D(:,:,j));
                    parfor e=layernumber(j)+1:layernumber(j+1)
                        stress(e,:,:,i)=squeeze(strain(e,:,:,i))*Dj;
                    end
                end
            end
        end
            
        function [ viscostress ] = viscoelasticstress3D(~,strain,...
                layernumber,Ei,T,D,v,t,n,elementnodes)
           for i = 1:size(D,3)  % layer number
                ee = exp(-t./T(i,:));  
                [cotese, cotes] = Results.cotes_calculation(Ei(i,:),T(i,:),t(i,:));
                ii = length(v) - size(D,3) + i; % layer index of vc layer
                cotesE = cotese/((1+v(ii))*(1-2*v(ii)))*...
                    [1-v(ii), v(ii), v(ii), 0, 0, 0;   
                     v(ii), 1-v(ii), v(ii), 0, 0, 0;
                     v(ii), v(ii), 1-v(ii), 0, 0, 0;
                     0, 0, 0, (1-2*v(ii))/2, 0, 0;
                     0, 0, 0, 0, (1-2*v(ii))/2, 0;
                     0, 0, 0, 0, 0, (1-2*v(ii))/2];
                Ej = zeros(length(T),1);
                for j = 1:length(T)
                    Ej(j) = Ei(i,j)*(ee(j)-1);
                end
            Cotesf = zeros(6,6,length(T));
            for j=1:length(T)
                Cotesf(:,:,j) = Ei(i,j)*cotes(j)/((1+v(ii))*(1-2*v(ii)))*...
                    [1-v(ii), v(ii), v(ii), 0, 0, 0;   
                     v(ii), 1-v(ii), v(ii), 0, 0, 0;
                     v(ii), v(ii), 1-v(ii), 0, 0, 0;
                     0, 0, 0, (1-2*v(ii))/2, 0, 0;
                     0, 0, 0, 0, (1-2*v(ii))/2, 0;
                     0, 0, 0, 0, 0, (1-2*v(ii))/2];
            end
            viscostress = zeros(size(elementnodes,1),size(elementnodes,2),6,n);
            Pn = zeros(size(elementnodes,1),size(elementnodes,2),6,n-1,length(T));
            j = 1; % initial state
            parfor e= layernumber(ii)+1:layernumber(ii+1)
                viscostress(e,:,:,j) = squeeze(strain(e,:,:,j))*...
                    squeeze(D(:,:,i))';
            end
            j = 2;
            temp1 = squeeze(viscostress(:,:,:,j-1));
            temp2 = squeeze((strain(:,:,:,j)-strain(:,:,:,j-1)));
            for e = layernumber(ii)+1:layernumber(ii+1)
                viscostress(e,:,:,j) = squeeze(temp1(e,:,:)) +...
                    squeeze(temp2(e,:,:))*cotesE;
                for jj = 1:length(T)
                    Pn(e,:,:,j-1,jj) = squeeze(temp2(e,:,:))*...
                        squeeze(Cotesf(:,:,jj));
                end
            end
            for j = 3:n % iteration
                temp1=squeeze(viscostress(:,:,:,j-1));
                temp2=squeeze((strain(:,:,:,j)-strain(:,:,:,j-1)));
                for e = layernumber(ii)+1:layernumber(ii+1)
                    temp3 = zeros(4,6);
                    for jj = 1:length(T)
                        temp3 = temp3+squeeze(Ej(jj)*Pn(e,:,:,j-2,jj));
                    end
                    viscostress(e,:,:,j) = squeeze(temp1(e,:,:)) +...
                        squeeze(temp2(e,:,:))*cotesE + temp3;
                    for jj = 1:length(T)
                        Pn(e,:,:,j-1,jj) = squeeze(ee(jj)*Pn(e,:,:,j-2,jj))+...
                            squeeze(temp2(e,:,:))*squeeze(Cotesf(:,:,jj));
                    end
                end
            end        
           end
        end
    end

    methods(Static)

        function [weights,locations] = GaussQuadrature()
        % Gauss quadrature for Q4 elements
               locations = ...
                  [ -0.57735 -0.57735;
                     0.57735 -0.57735;
                     0.57735  0.57735;
                    -0.57735  0.57735];
                weights = [1;1;1;1];    
        end

        function [JacobianMatrix,XYDerivatives] = ...
                Jacobian(nodeCoordinates,naturalDerivatives)   
            JacobianMatrix = nodeCoordinates'*naturalDerivatives;                   
            XYDerivatives = naturalDerivatives/JacobianMatrix;
        end

        function [shape,naturalderivatives] = ShapeFunctionQ4(xi,eta)
        shape = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...
                     (1+xi)*(1+eta);(1-xi)*(1+eta)];
        naturalderivatives = 1/4*[-(1-eta),-(1-xi);
                                  (1-eta),-(1+xi);
                                  (1+eta),(1+xi);
                                  -(1+eta),(1-xi)];         
        end
    
        function [cotes, temp] = cotes_calculation(K,T,t)
            temp=zeros(1,length(T));
            for i=1:length(temp)
                temp(i)=1/90*(7*exp(-t/T(i))+...
                    32*exp(-0.75*t/T(i))+12*exp(-0.5*t/T(i))...
                    +32*exp(-0.25*t/T(i))+7);
            end
            cotes=K(end)+sum(K(1:end-1).*temp);
        end
    end
end

