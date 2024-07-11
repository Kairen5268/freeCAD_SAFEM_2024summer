classdef Contact_stress
    %CONTACT_STRESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pressure
        len
    end
    
    methods
        function obj = Contact_stress(P,S,tire,process)
            % generate the contact stress of each rib
            % Ref: Analytical approach for predicting three-dimensional tire-pavement contact load
            % only considering vertical stress
            % Transportation Research Record, 2014, Hernandez, Jaime A.et. al
            switch tire.tire_name
                case 'WBT 445'
                    obj.pressure = zeros(length(tire.rib_width),201,process.step);
                    obj.len = zeros(length(tire.rib_width),201,process.step);
                    for i = 1:process.step
                        [obj.pressure(:,:,i), obj.len(:,:,i)]...
                            = obj.WBT(P*process.loading_factor(i),S,tire);
                    end

                otherwise

            end
        end
        
        
          
     
    end

    methods (Static)

      function [pressure,len] = WBT(P,S,tire)
        % 445 wide-base tire with 8 ribs% input tire related parameters
        d1=[76.2 117.7 140.9 147.0 147.9 143.7 121.9 77.9];
        d2=[3.573 2.182 2.039 2.090 2.148 1.972 1.872 3.144];
        k=[1.229 1.093 1.045 1.526 1.461 0.901 1.156 1.050;
            7.824 8.159 10.044 10.748 8.602 8.698 9.092 10.102;
            -0.674 -0.567 -0.608 -1.152 -1.009 -0.476 -0.726 -0.745];
        c=[0.042 0.056 0.063 0.096 0.080 0.116 0.025 -0.008;
            1.684 -0.366 -0.754 -0.502 -0.729 -0.767 -0.845 -0.009;
            -0.073 -0.025 0.025 -0.090 0.003 -0.158 0.108 0.097;
            -8.978 1.280 3.945 2.037 3.243 4.343 5.290 5.313;
            0.016 0.038 -0.017 0.081 0.006 0.127 -0.063 -0.075];
        % contact length
        l=d1+P*d2;
        % parameters for n and alpha
        n=k(1,:)+P*1e-3*k(2,:)+S*1e-3*k(3,:);
        alpha=c(1,:)+P*1e-3*c(2,:)+S*1e-3*c(3,:)+(P*1e-3)^2*c(4,:)+(S*1e-3)^2*c(5,:);
        % contact stress
        q=zeros(8,101);
        i=1;
        for x=0:0.01:1
            for j=1:8
                q(j,i)=alpha(j)*(P/l(j))*(1+1/(2*n(j)))*(1-x^(2*n(j)));
            end
            i=i+1;
        end
        % check
        loading=zeros(1,8);
        for i=1:8
            loading(i)=sum(l(i)/100*q(i,:));
        end
        temp=sum(loading);
        factor=P/temp;
        q_modified=q*factor;
        len=zeros(8,201);
        pressure=zeros(8,201);
        width=tire.rib_width;
        for i=1:8
            len(i,:)=-l(i)/2:l(i)/200:l(i)/2;
            % unit MPa
            pressure(i,:)=[flip(q_modified(i,2:end)) q_modified(i,:)]/width(i);
        end


      end

    end
end

