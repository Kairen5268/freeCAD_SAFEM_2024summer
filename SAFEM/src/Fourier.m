classdef Fourier
    %FOURIER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        loadfourier3D
        factor_corr
        Lz
        L
    end
    
    methods
        function obj = Fourier(Lz,tire,process,contact_stress,inputArg4)
            %FOURIER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Lz = Lz;
            if nargin >= 5  
                obj.factor_corr = inputArg4;  
            else
                obj.factor_corr = 0.95; 
            end            
            obj.loadfourier3D = obj.Loading3Dfourier(tire,...
                process,contact_stress); 
            obj.L = unique(obj.loadfourier3D(:,1,:));
        end
        
        function loadfourier3D = Loading3Dfourier(obj,...
                tire,process,contact_stress)
            % fourier dimension and corresponding weight
            rib_num = length(tire.rib_width);
            loadfourier3D = zeros(1001,6*process.step+1,rib_num);
            templ=zeros(1,rib_num);
            for i = 1:rib_num
                temp = Fourier_rib(obj,process,...
                    contact_stress,i);
                loadfourier3D(1:size(temp,1),:,i) = temp;
                templ(i)=size(temp,1);
            end
            loadfourier3D(max(templ)+1:end,:,:)=[];
        end

        function loadfourier3Drib = Fourier_rib(obj,process,contact_stress,rib_num)
            % calculate fourier for each rib loading
            lz = obj.Lz;
            factor = obj.factor_corr;
            lenz = squeeze(contact_stress.len(rib_num,:,:))/1000; % m
            pressurez = squeeze(contact_stress.pressure(rib_num,:,:));
            zcenter = process.z_coordinate;
            fourier_rib = zeros(1001,2*process.step+1);
            fourier_rib(:,1) = 0:1000;
            temp1=zeros((2*lz/0.01+1)+1,process.step);
            temp2=zeros((2*lz/0.01+1)+1,process.step);
            parfor i = 1:process.step
                load_raw=[0 0; lenz(:,i)+zcenter(i) pressurez(:,i);  lz 0];
                load_data = Fourier.temp_interp_loading(load_raw,lz);
                [temp1(:,i),~,temp2(:,i)] = Fourier.temp_Loading_FFT(load_data,factor);
            end
            N=length(-lz:0.01:lz);
            for i=1:process.step
                entire_index=temp1(1:nnz(temp1(:,i)),i);
                F_entire=temp2(:,i);
                for j=1:length(entire_index)
                    if (entire_index(j)-(N+1)/2-1)>0 % sin
                       L=(entire_index(j)-(N+1)/2-1);
                       fourier_rib(L+1,2*i+1)=F_entire(entire_index(j));
                    else                             % cos
                       L=(entire_index(j)-1);
                       fourier_rib(L+1,2*i)=F_entire(entire_index(j));
                    end
                end
            end
            fourier_rib(all(fourier_rib(:,2:end)==0,2),:)=[];
            fourier_rib(:,2:end)=2/N*fourier_rib(:,2:end);
            fourier_rib(1,2:end)=fourier_rib(1,2:end)/2;
            % 3D loading fourier
            loadfourier3Drib=zeros(1001,6*process.step+1);
            loadfourier3Drib(:,1)=0:1000;
            loadfourier3Drib(fourier_rib(:,1)+1,4:6:end)=fourier_rib(:,2:2:end);
            loadfourier3Drib(fourier_rib(:,1)+1,5:6:end)=fourier_rib(:,3:2:end);
            loadfourier3Drib(all(loadfourier3Drib(:,2:end)==0,2),:)=[];
            end
        end

    methods (Static)

        function [load_data] = temp_interp_loading(load_raw,Lz)
            %% generate loading interpolation (0.01 m)
            xq=0:0.01:Lz;
            vq=interp1(load_raw(:,1),load_raw(:,2),xq);
            % sin function
            temp1=[flip(-xq(2:end)) xq]';
            temp2=[flip(-vq(2:end)) vq]';
            load_data=[temp1 temp2];
        end

        function [entire_index,N,F_entire] = temp_Loading_FFT(load_data,factor_corr)
            %% 2021.3.28 nanjing, SAPAVE 2.0 (finished)
            % obtain the weighting Fourier series term of the loading function
            % superposition of cos and sin, verification by cross-correlation
            % input: load_data: load (coordinate z vs. load) from excel
            %        factor_corr (cross-correlation factor)
            %% discrete loading data
            z=load_data(:,1); load_original=load_data(:,2);
            % obtain the load function for the FFT (FFT will generate inversion !!!)
            load=load_original;
            load(2:length(load)-1)=flip(load_original(2:length(load)-1));
            %% fast fourier transform
            % keep in mind: the FFT generates the mirror image of Fourier coefficient
            % fftshift function should be conducted to obtain the true coefficient
            N=length(z); fourier=fftshift(fft(load));
            F_real=real(fourier((N+1)/2:N)); F_imag=imag(fourier((N+1)/2:N));
            %% coefficient weighting method for inverse transform
            cos_temp=abs(F_real'); sin_temp=abs(F_imag'); F_entire=[F_real' F_imag'];
            entire_temp=[cos_temp sin_temp]; 
            entire_index=zeros(length(entire_temp),1); entire_value=zeros(length(entire_temp),1);
            load_weight=zeros(N,1);
            temp=zeros(N+1,1);
            % check by cross correlation, note the NaN when vector is zero
            i=0;
            while (xcorr(load_original,load_weight,0,'normalized')<factor_corr-1e-5||...
                    isnan(xcorr(load_original,load_weight,0,'normalized')))
                i=i+1;
                [m,index]=max(entire_temp);
                entire_index(i)=index; entire_value(i)=m; 
                entire_temp(index)=0;
                for n=1:N
                    for k=1:i
                        if entire_index(k)<=(N+1)/2
                           temp(entire_index(k))=cos(2*pi*(entire_index(k)-1)*n/N);
                        else
                           temp(entire_index(k))=sin(2*pi*(entire_index(k)-(N+1)/2-1)*n/N);
                        end
                    end
                    load_weight(n)=(-F_real(1)+2*F_entire(entire_index(1:i))*temp(entire_index(1:i)))/N; % *2!!!
                end
            end
            %% maxmum loading modification
            factor=max(load_original)/max(load_weight);
            F_entire=F_entire*factor;           
        end

    end
  
end

