classdef Solution
    %SOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Displacement
        Velocity
        Acceleration
        depthtemp
        visco_nodeindex
    end
    
    methods
        function obj = Solution(dimension,tire,meshing,process,contactstress,...
                fourier,material,boundary,interface,temperature)
            %SOLUTION Construct an instance of this class

            L=unique(fourier.loadfourier3D(:,1,:));
            Displacement=zeros(meshing.GDof,process.step,length(L));
            Velocity=zeros(meshing.GDof,process.step,length(L));
            Acceleration=zeros(meshing.GDof,process.step,length(L));
            % interation
            GDof = meshing.GDof;
            numberelements = length(meshing.element_node);
            elementdata = meshing.element_node;
            numbernode = meshing.GDof/3;
            nodecoordinate = meshing.nodecoordinate;
            layernumber = [0 meshing.layer_element];
            Loading_index = meshing.loading_index;
            D = material.D;
            rho = material.rho;
            ray = material.rayleigh;
            k1 = material.k1;
            k2 = material.k2;
            Ki = material.Ki;
            Gi = material.Gi;
            T = material.T;
            WLC = material.WLF;
            n = process.step;
            t = process.time;
            Lz = fourier.Lz;
            loadfourier3D = fourier.loadfourier3D;
            Bottom = boundary.Bottom;
            Vertical = boundary.Vertical;
            interfacestiff = interface.stiffness;
            temperaturedepth = temperature.one_dimension;
            tempindex = layernumber(end-size(Gi,1)); %*
            visco_nodeindex = elementdata(tempindex+1,1); %*
            obj.depthtemp = -nodecoordinate(visco_nodeindex:end,2)+...
            nodecoordinate(end,2); %*
            depthtemp = obj.depthtemp; %*
            tempindex = layernumber(end-size(Gi,1));
            obj.visco_nodeindex = elementdata(tempindex+1,1);
            parfor i = 1:length(L)
                l = L(i);
                [Displacement(:,:,i),Velocity(:,:,i),Acceleration(:,:,i)] = ...
                    Solution.Fouriersolution(GDof,numberelements,elementdata,numbernode,...
                    nodecoordinate,D,Lz,l,layernumber,rho,ray,n,k1,k2,Ki,Gi,T,...
                    loadfourier3D,t,Loading_index,Bottom,Vertical,...
                    interfacestiff,temperaturedepth,WLC,depthtemp);
            end
            obj.Displacement = Displacement;
            obj.Velocity = Velocity;
            obj.Acceleration = Acceleration;
        end
    end

    methods(Static)
        %% 
        function [Displacement,Velocity,Acceleration] = Fouriersolution(...
            GDof,numberelements,elementdata,numbernodes,...
            nodecoordinates,D,Lz,l,layernumber,rho,ray,...
            n,k1,k2,Ki,Gi,Ti,loadfourier3D,t,...
            Loading_index,Bottom,Vertical,...
            interface,temperature_data,WLC,depthtemp)
        % preallocation variables 
            Displacement=zeros(GDof,n);
            Velocity=zeros(GDof,n);
            Acceleration=zeros(GDof,n);
        % viscoelastic layer node index
        tempindex = layernumber(end-size(Gi,1));
        visco_nodeindex = elementdata(tempindex+1,1);
        % material constant parameters
        % elastic martix
        [Stiffness,StiffnessD4,Mass,Damping] = Solution.formmartixfourier(...
            GDof,numberelements,...
            elementdata,numbernodes,nodecoordinates,...
            D,Lz,l,layernumber,rho,ray);
        % Cotes (Gi+Ki) martix
        Temp1=Solution.formCotesKfourier(GDof,elementdata,numbernodes,...
            nodecoordinates,k1,k2,Ki,Gi,t,Ti,Lz,l,...
            layernumber,temperature_data,WLC);
        % loading function, (GDof,2*(n+1)) cos & sin
        Loading=Solution.temp_loadgeneration(GDof,n,l,loadfourier3D,...
            Loading_index,Lz);
        Loading=[Loading(:,1)*(0:0.1:1) Loading(:,2:end)];
        % original state
        % viscoelastic materials is regared as elastic
        Fload=squeeze(Loading(:,1));
        Ktemp=Stiffness+StiffnessD4+Bottom+interface;
        Ftemp=Fload;
        Displacement(:,1)=Solution.solution_CPU(Ktemp,Ftemp);
        % dynamic viscoelastic analysis (virtual work principle)
        % viscoelastic stress temp
        Fviscotemp=zeros(GDof,n);
        Fviscotemp(:,1)=StiffnessD4*squeeze(Displacement(:,1));
        % first step
        i=1;
        % mass force, temp for parfor cycle
            Fmass=Mass*squeeze(4/(t^2)*Displacement(:,i)+4/t*Velocity(:,i)+Acceleration(:,i));
            Fdamp=Damping*squeeze(2/t*Displacement(:,i)+Velocity(:,i));
            Fload=squeeze(Loading(:,i+1)); % -
            Fvertical=Vertical*squeeze(2/t*Displacement(:,i)+Velocity(:,i));
            F0visco=Fviscotemp(:,1);
            Fvisco=Temp1*squeeze(Displacement(:,i));
            Ktemp=Mass*4/t^2+Damping*2/t+Stiffness+Temp1+Bottom+Vertical*2/t+interface;
            Ftemp=Fload+Fvertical+Fmass+Fdamp-F0visco+Fvisco;
            Displacement(:,i+1)=Solution.solution_CPU(Ktemp,Ftemp);
            Velocity(:,i+1)=2*(Displacement(:,i+1)-Displacement(:,i))/t-Velocity(:,i);
            Acceleration(:,i+1)=4*(Displacement(:,i+1)-Displacement(:,i))/t^2-...
                4*Velocity(:,i)/t-Acceleration(:,i);
            Fviscotemp(:,i+1)=F0visco+Temp1*squeeze(Displacement(:,i+1)-Displacement(:,i));
        % interation cycle
        % Tempcotes Pn+Qn
        Pnforce=zeros(GDof,length(Ti(:)));
        displacementtemp=squeeze(Displacement(:,2)-Displacement(:,1));
        for i=1:length(Ti(:))
            [px,py]=find(Ti==Ti(i)); % Pn location
            temp=size(Ti,1)-px+1;
            Temp=Solution.formtempcotesfourier(GDof,elementdata,...
                numbernodes,nodecoordinates,t,Ti(i),Gi(px,py),Ki(px,py),k1,k2,Lz,l,...
                layernumber(end-temp:end-temp+1),temperature_data,WLC);
            eval([['TempcotesPQ_',num2str(i)],'= Temp;']);
                Pnforce(:,i)=eval(['TempcotesPQ_',num2str(i)])*displacementtemp;
        end
        for i=2:n-1         
            Fmass=Mass*squeeze(4/(t^2)*Displacement(:,i)+4/t*Velocity(:,i)+Acceleration(:,i));
            Fdamp=Damping*squeeze(2/t*Displacement(:,i)+Velocity(:,i));
            Fload=squeeze(Loading(:,i+1)); % -
            Fvertical=Vertical*squeeze(2/t*Displacement(:,i)+Velocity(:,i));
            F0visco=Fviscotemp(:,i);
            Fvisco=Temp1*squeeze(Displacement(:,i))-sum(Pnforce,2); % 右边
            Ktemp=Mass*4/t^2+Damping*2/t+Stiffness+Temp1+Bottom+Vertical*2/t+interface;
            Ftemp=Fload+Fvertical+Fmass+Fdamp-F0visco+Fvisco;
            Displacement(:,i+1)=Solution.solution_CPU(Ktemp,Ftemp);
            Velocity(:,i+1)=2*(Displacement(:,i+1)-Displacement(:,i))/t-Velocity(:,i);
            Acceleration(:,i+1)=4*(Displacement(:,i+1)-Displacement(:,i))/t^2-4*Velocity(:,i)/t-Acceleration(:,i);
            Fviscotemp(:,i+1)=F0visco+Temp1*squeeze(Displacement(:,i+1)-Displacement(:,i))+sum(Pnforce,2);
            displacementtemp=squeeze(Displacement(:,i+1)-Displacement(:,i));
            for j=1:size(Pnforce,2)
                   Pnforce(:,j)=Solution.temp_Pnforceinteration(numbernodes,nodecoordinates,...
                    t,Ti(j),Pnforce(:,j),eval(['TempcotesPQ_',num2str(j)]),...
                    displacementtemp,temperature_data,WLC,visco_nodeindex,elementdata,depthtemp);
            end
        end
        end
        
        function [stiffness_e,stiffness_v,mass,damping] = formmartixfourier(...
                Gdof,~,elementdata,...
            numbernodes,nodecoordinates,D,Lz,l,layernumber,rho,ray)
        % three martix
        damping=zeros(Gdof);
        stiffness=zeros(Gdof);
        mass=zeros(Gdof);
        [gaussweights,gausslocations]=Solution.GaussQuadrature();
        gammal=l*pi/Lz;
        for i=1:length(rho)  % layer number
        for e=(layernumber(i)+1):layernumber(i+1)
            indice=elementdata(e,:);
            elementdof=[indice indice+numbernodes indice+2*numbernodes];
            ndof=length(indice);
        for q=1:size(gaussweights,1) % Gauss point
                gausspoint=gausslocations(q,:);
                xi=gausspoint(1);eta=gausspoint(2);
                [shapefunction,naturalderivatives]=...
                    Solution.ShapeFunctionQ4(xi,eta);
                [jacob,XYderivatives]=Solution.Jacobian(...
                    nodecoordinates(indice,:),naturalderivatives);
        % B matrix
        B=zeros(6,3*ndof);
        B(1,1:ndof)=XYderivatives(:,1)';
        B(2,ndof+1:2*ndof)=XYderivatives(:,2)';
        B(3,2*ndof+1:3*ndof)=-shapefunction*gammal;
        B(4,1:ndof)=XYderivatives(:,2)';
        B(4,ndof+1:2*ndof)=XYderivatives(:,1)';
        B(5,ndof+1:2*ndof)=shapefunction*gammal;
        B(5,2*ndof+1:3*ndof)=XYderivatives(:,2)';
        B(6,1:ndof)=shapefunction*gammal;
        B(6,2*ndof+1:3*ndof)=XYderivatives(:,1)';
        % matrix calculation
        stiffness(elementdof,elementdof)=stiffness(elementdof,elementdof)+B'*D(:,:,i)*B*gaussweights(q)*det(jacob)*Lz/2;
        mass(indice,indice)=mass(indice,indice)+shapefunction*shapefunction'*gaussweights(q)*rho(i)*det(jacob)*Lz/2;
        mass(indice+numbernodes,indice+numbernodes)=mass(indice+numbernodes,indice+numbernodes)+...
        shapefunction*shapefunction'*gaussweights(q)*rho(i)*det(jacob)*Lz/2;
        mass(indice+2*numbernodes,indice+2*numbernodes)=mass(indice+2*numbernodes,indice+2*numbernodes)+...
        shapefunction*shapefunction'*gaussweights(q)*rho(i)*det(jacob)*Lz/2;
        end
        end
        if i==length(ray)
        damping=ray(1,1)*mass+ray(1,2)*stiffness;
        stiffness_e=stiffness;
        end
        end
        stiffness_e=sparse(stiffness_e);
        stiffness_v=sparse(stiffness-stiffness_e);
        mass=sparse(mass);
        damping=sparse(damping);
        end
        
        function [stiffness] = formCotesKfourier(...
                Gdof,elementdata,numbernodes,...
                nodecoordinates,k1,k2,K,G,t,TT,Lz,l,...
                layernumber,temperature_data,WLC)
        % form Cotes(Ki+Gi) matrix, two viscoelastic layers
        stiffness=zeros(Gdof);
        [gaussweights,gausslocations] = Solution.GaussQuadrature();
        gammal=l*pi/Lz;
        for i=1:size(K,1) % layer number
            Ki=K(i,:); 
            Gi=G(i,:); 
            Ti=TT(i,:);
        % element number
        for e=(layernumber(end-size(K,1)+i-1)+1):layernumber(end-size(K,1)+i)
            indice=elementdata(e,:);
            % temperature parameters (element center as temperature reference)
            depth=nodecoordinates(end,2)-mean(nodecoordinates(indice,2));
            [Tshift]=Solution.temperature(temperature_data,depth,WLC);
            T=Ti*Tshift;
            % Ki+Gi
            g=2*Solution.cotes_calculation(Gi,T,t);         
            k=3*Solution.cotes_calculation(Ki,T,t);   
            D=g*k1+k*k2;          
            elementdof=[indice indice+numbernodes indice+2*numbernodes];
            ndof=length(indice);
            % Gauss point iteration
            for q=1:size(gaussweights,1)
                gausspoint=gausslocations(q,:);
                xi=gausspoint(1);eta=gausspoint(2);
                [shapefunction,naturalderivatives]=Solution.ShapeFunctionQ4(xi,eta);
                [jacob,XYderivatives]=Solution.Jacobian(nodecoordinates(indice,:),naturalderivatives);
                % B matrix
                B=zeros(6,3*ndof);
                B(1,1:ndof)=XYderivatives(:,1)';
                B(2,ndof+1:2*ndof)=XYderivatives(:,2)';
                B(3,2*ndof+1:3*ndof)=-shapefunction*gammal;
                B(4,1:ndof)=XYderivatives(:,2)';
                B(4,ndof+1:2*ndof)=XYderivatives(:,1)';
                B(5,ndof+1:2*ndof)=shapefunction*gammal;
                B(5,2*ndof+1:3*ndof)=XYderivatives(:,2)';
                B(6,1:ndof)=shapefunction*gammal;
                B(6,2*ndof+1:3*ndof)=XYderivatives(:,1)';
                stiffness(elementdof,elementdof)=stiffness(elementdof,elementdof)+...
                    B'*D*B*gaussweights(q)*det(jacob)*Lz/2;
            end
        end
        end
        stiffness=sparse(stiffness);
        end

        function [stiffness] = formtempcotesfourier(Gdof,elementdata,...
        numbernodes,nodecoordinates,t,Ti,Gi,Ki,k1,k2,Lz,l,...
        layernumber,temperature_data,WLC)
        % Pn+Qn matrix, considering Ki
        stiffness=zeros(Gdof);
        [gaussweights,gausslocations]=Solution.GaussQuadrature();
        gammal=l*pi/Lz;
        % element number
        for e=(layernumber(1)+1):layernumber(2)
            indice=elementdata(e,:);
            % temperature parameters
            depth=nodecoordinates(end,2)-mean(nodecoordinates(indice,2));
            [Tshift]=Solution.temperature(temperature_data,depth,WLC);
            T=Ti*Tshift; 
            cotes=1/90*(7*exp(-t/T)+32*exp(-0.75*t/T)+12*exp(-0.5*t/T)+32*exp(-0.25*t/T)+7);
            D=(exp(-t/T)-1)*(2*Gi*k1+3*Ki*k2)*cotes;        
            elementdof=[indice indice+numbernodes indice+2*numbernodes];
            ndof=length(indice);
            % Gauss point
            for q=1:size(gaussweights,1)
                gausspoint=gausslocations(q,:);
                xi=gausspoint(1);eta=gausspoint(2);
                [shapefunction,naturalderivatives]=Solution.ShapeFunctionQ4(xi,eta);
                [jacob,XYderivatives]=Solution.Jacobian(nodecoordinates(indice,:),naturalderivatives);
                % B matrix
                B=zeros(6,3*ndof);
                B(1,1:ndof)=XYderivatives(:,1)';
                B(2,ndof+1:2*ndof)=XYderivatives(:,2)';
                B(3,2*ndof+1:3*ndof)=-shapefunction*gammal;
                B(4,1:ndof)=XYderivatives(:,2)';
                B(4,ndof+1:2*ndof)=XYderivatives(:,1)';
                B(5,ndof+1:2*ndof)=shapefunction*gammal;
                B(5,2*ndof+1:3*ndof)=XYderivatives(:,2)';
                B(6,1:ndof)=shapefunction*gammal;
                B(6,2*ndof+1:3*ndof)=XYderivatives(:,1)';
                stiffness(elementdof,elementdof)=stiffness(elementdof,elementdof)+...
                    B'*D*B*gaussweights(q)*det(jacob)*Lz/2;
            end
        end
        stiffness=sparse(stiffness);
        end

        function [cotes] = cotes_calculation(K,T,t)
        temp=zeros(1,length(T));
        for i=1:length(temp)
            temp(i)=1/90*(7*exp(-t/T(i))+32*exp(-0.75*t/T(i))+12*exp(-0.5*t/T(i))...
                +32*exp(-0.25*t/T(i))+7);
        end
        cotes=K(end)+sum(K(1:end-1).*temp);
        end
        
        function [Pnforce] = temp_Pnforceinteration(...
            numbernodes,nodecoordinates,t,Ti,Pnforce,K,...
            displacement,temperature_data,WLC,visco_nodeindex,elementdata,depthtemp)
        % calculation Pnforce interation
        Pntemp=zeros(length(Pnforce),1);
        depth=mean(depthtemp);
        
        for i=visco_nodeindex:numbernodes
            % temp=find(elementdata(:,4)==i);
            % indice=elementdata(e,:);
            % % temperature parameters (element center as temperature reference)
            % depth=nodecoordinates(end,2)-mean(nodecoordinates(indice,2));

           [Tshift]=Solution.temperature(temperature_data,depth,WLC);
           T = Ti*Tshift; 
           
           Pntemp(i+[0 1 2]*numbernodes)=...
               Pnforce(i+[0 1 2]*numbernodes)*exp(-t/T);
        end
        Pnforce=Pntemp+K*displacement;
        
        end

        function [Tshift] = temperature(temperature_data,depth,WLC)
        % 2022.1.14 Piscataway, SAPAVE 2.0 (finished)
        % obtain the temperature shift factor for each element
        % Input: element temperature T & reference temperature Tref
        %        equation parameter (vector, containing all required parameters)
        % further study can consider different temperature model
        T=50;
        for i=1:length(temperature_data)
            if (abs(temperature_data(i,1)-depth*1000)<1)
               T=temperature_data(i,2);
               break
            end
        end
        % obtain the shift factor
        Tref=50;
        lg_Tshift=-WLC(1)*(T-Tref)/(WLC(2)+(T-Tref));
        Tshift=10^lg_Tshift;  
        end

        function [loading] = temp_loadgeneration(GDof,n,l,...
            loadfourier3D,Loading_index,Lz)
        % generate loading function considering dynamic loads
        loading=zeros(GDof,n); % sin
        % vertical loading
        for i=1:length(Loading_index)/3 % rib number
            if l==0
                loading(Loading_index(3*(i-1)+1:3*i,1)+GDof/3,:)=Lz/2*squeeze(Loading_index(3*(i-1)+1:3*i,2))*...
                    ((-1)^l*squeeze(loadfourier3D(1,5:6:end,i))); %sin
            else 
                j=find(loadfourier3D(:,1,i)==l);
                if isempty(j)
                else
                    loading(Loading_index(3*(i-1)+1:3*i,1)+GDof/3,:)=Lz/2*squeeze(Loading_index(3*(i-1)+1:3*i,2))*...
                    ((-1)^l*squeeze(loadfourier3D(j,5:6:end,i)));
                end
            end
        end
        loading=-1e6*loading; % change direction and MPa to Pa
        end

        function [displacements] = solution_CPU(stiffness,force)        
        displacements=stiffness\force;
        end

        function [weights,locations]=GaussQuadrature()
        % Gauss quadrature for Q4 elements
               locations=...
                  [ -0.57735 -0.57735;
                     0.57735 -0.57735;
                     0.57735  0.57735;
                    -0.57735  0.57735];
                weights=[1;1;1;1];    
        end
        
        function [shape,naturalderivatives]=ShapeFunctionQ4(xi,eta)
        shape=1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
        naturalderivatives=1/4*[-(1-eta),-(1-xi);
                                (1-eta),-(1+xi);
                                (1+eta),(1+xi);
                                -(1+eta),(1-xi)] ;         
        end
        
        function [JacobianMatrix,XYDerivatives]=Jacobian(...
                nodeCoordinates,naturalDerivatives)   
            JacobianMatrix=nodeCoordinates'*naturalDerivatives;                   
            XYDerivatives=naturalDerivatives/JacobianMatrix;
        end

    end
end

