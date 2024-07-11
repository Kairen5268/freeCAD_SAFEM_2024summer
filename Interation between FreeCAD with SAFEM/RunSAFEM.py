import matlab.engine
import numpy as np

# Start MATLAB engine
eng = matlab.engine.start_matlab()

# Function to convert Python lists to MATLAB arrays
def to_matlab_array(py_list):
    return matlab.double(np.array(py_list).tolist())

# Dimensions and material properties from FreeCAD
dimensions = [2, 0.5, 0.2, 0.1]
material_properties = {
    "rho": [2100, 2200, 2300, 2043],
    "v": [0.45, 0.25, 0.25, 0.35],
    "E": [103e6, 517e6, 1723e6],
    "ray": [[1.04, 5.59e-3], [1.04, 5.59e-3], [1.04, 5.59e-3]],
    "Ei": [6709, 3531, 4225, 3217, 2078, 944, 387, 409],
    "T": [1.001e-3, 1.001e-2, 1.001e-1, 1.001, 1.001e1, 1.001e2, 1.001e3],
    "WLF": [50, 315]
}

# Convert the lists to MATLAB arrays
dimensions_m = to_matlab_array(dimensions)
rho_m = to_matlab_array(material_properties["rho"])
v_m = to_matlab_array(material_properties["v"])
E_m = to_matlab_array(material_properties["E"])
ray_m = to_matlab_array(material_properties["ray"])
Ei_m = to_matlab_array(material_properties["Ei"])
T_m = to_matlab_array(material_properties["T"])
WLF_m = to_matlab_array(material_properties["WLF"])

# Run the MATLAB script with the data
eng.workspace['dimension'] = dimensions_m
eng.workspace['rho'] = rho_m
eng.workspace['v'] = v_m
eng.workspace['E'] = E_m
eng.workspace['ray'] = ray_m
eng.workspace['Ei'] = Ei_m
eng.workspace['T'] = T_m
eng.workspace['WLF'] = WLF_m

eng.eval("dimension = Dimension(dimension,5);", nargout=0)
eng.eval("tire = Tire('WBT 445',1);", nargout=0)
eng.eval("meshing = Meshing(dimension,tire);", nargout=0)
eng.eval("process = Process(20,20,0.005,1);", nargout=0)
eng.eval("contactstress = Contact_stress(47,724,tire,process);", nargout=0)
eng.eval("fourier = Fourier(10,tire,process,contactstress,0.95);", nargout=0)
eng.eval("material = Material(rho,v,E,ray,Ei,T,WLF);", nargout=0)
eng.eval("boundary = Boundary(material,meshing,fourier);", nargout=0)
eng.eval("interface = Interface(meshing,fourier);", nargout=0)
eng.eval("temperature = Temperature('input',[12.7 36.1; 38.1 55.2; 76.2 54.1; 127 46.9]);", nargout=0)
eng.eval("solution = Solution(dimension,tire,meshing,process,contactstress, fourier,material,boundary,interface,temperature);", nargout=0)
eng.eval("results = Results('fixed_location',2,meshing,process,fourier,solution,material);", nargout=0)
eng.eval("ResultsFigure('Cross-section longitudinal strain',results,meshing);", nargout=0)

# Close MATLAB engine
eng.quit()
