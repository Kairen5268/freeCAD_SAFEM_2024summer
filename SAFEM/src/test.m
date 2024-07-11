%% 07112024 test for FreeCAD 2024 summer
dimension = Dimension([2 0.5 0.2 0.1],5);
tire = Tire('WBT 445',1);
meshing = Meshing(dimension,tire);
% Figure.mesh(meshing);
process = Process(20,20,0.005,1);
contactstress = Contact_stress(47,724,tire,process);
% Figure.contact_stress(contactstress,tire,10);
fourier = Fourier(10,tire,process,contactstress,0.95);
% Figure.loading_fourier(fourier,contactstress,process,4,2);
rho=[2100 2200 2300 2043];
v=[0.45 0.25 0.25 0.35];
E=[103e6 517e6 1723e6];
ray=[1.04 5.59e-3; 1.04 5.59e-3; 1.04 5.59e-3];
Ei=[6709 3531 4225 3217 2078 944 387 409]*1e6;   % wearing course     
T=[1.001e-3 1.001e-2 1.001e-1 1.001 1.001e1 1.001e2 1.001e3];
WLF = [50 315];
material = Material(rho,v,E,ray,Ei,T,WLF);
boundary = Boundary(material,meshing,fourier);
interface = Interface(meshing,fourier);
temperature = Temperature('input',[12.7 36.1; 38.1 55.2; 76.2 54.1; 127 46.9]);
% Figure.temperature_1D(temperature,dimension,material);
solution = Solution(dimension,tire,meshing,process,contactstress,...
                fourier,material,boundary,interface,temperature);
%%
results = Results('fixed_location',2,meshing,process,fourier,solution,material);
%%
ResultsFigure('Cross-section longitudinal strain',results,meshing);                         
      