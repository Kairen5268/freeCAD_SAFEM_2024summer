classdef Material
    %MATERIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rho
        v
        E
        D
        rayleigh
        Ei
        Ki
        Gi
        T
        k1
        k2
        WLF
    end
    
    methods
        function obj = Material(rho,v,E,ray,Ei,T,WLF)
            obj.rho = rho;
            obj.v = v;
            obj.E = E;
            obj.rayleigh = ray;
            obj.Ei = Ei;
            obj.T = T;
            obj.D = zeros(6,6,length(v));
            obj.WLF = WLF;
            for i=1:length(E)
                obj.D(:,:,i)=E(i)/((1+v(i))*(1-2*v(i)))*...
                    [1-v(i ), v(i) , v(i) ,  0   ,    0   ,   0    ;   
                     v(i) ,1-v(i) , v(i) ,  0    ,    0   ,   0    ;  
                     v(i) ,  v(i) ,1-v(i),  0    ,    0   ,   0    ; 
                     0  ,  0 , 0 ,(1-2*v(i))/2,    0    ,    0     ;
                     0  ,  0 , 0 ,  0      ,(1-2*v(i))/2,    0     ;
                     0  ,  0 , 0 ,  0      ,    0    ,(1-2*v(i))/2]; 
            end
            obj.Ki=zeros(size(Ei)); obj.Gi = obj.Ki;
            % instant mudulus
            for i=length(E)+1:length(v)
                obj.Ki(i-length(E),:)=Ei(i-length(E),:)/(3*(1-2*v(i)));
                obj.Gi(i-length(E),:)=Ei(i-length(E),:)/(2*(1+v(i)));
                obj.D(:,:,i)=sum(Ei(i-length(E),:))/((1+v(i))*(1-2*v(i)))*...
                    [1-v(i), v(i) , v(i) ,  0      ,    0    ,    0     ;   
                     v(i) ,1-v(i) , v(i) ,  0      ,    0    ,    0     ;  
                     v(i) ,  v(i) ,1-v(i),  0      ,    0    ,    0     ; 
                     0  ,  0 , 0 ,(1-2*v(i))/2,    0    ,    0     ;
                     0  ,  0 , 0 ,  0      ,(1-2*v(i))/2,    0     ;
                     0  ,  0 , 0 ,  0      ,    0    ,(1-2*v(i))/2];
            end
            % generate viscoelastic material parameters
            obj.k1 = [2/3 -1/3 -1/3 0 0 0;
                      -1/3 2/3 -1/3 0 0 0;
                      -1/3 -1/3 2/3 0 0 0;
                      0 0 0 1/2 0 0 ;
                      0 0 0 0 1/2 0 ;
                      0 0 0 0 0 1/2 ];% deviator matrix
            obj.k2 = [1/3 1/3 1/3 0 0 0;  
                    1/3 1/3 1/3 0 0 0;
                    1/3 1/3 1/3 0 0 0;
                    0 0 0 0 0 0;
                    0 0 0 0 0 0;
                    0 0 0 0 0 0];% principal matrix
        end
    end
end

