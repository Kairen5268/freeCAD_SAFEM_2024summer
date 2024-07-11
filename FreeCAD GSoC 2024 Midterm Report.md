## FreeCAD GSoC 2024 Midterm Report: Implementation of Semi-Analytical Finite Element in FreeCAD for Pavement Infrastructure Modeling

## Kairen Shen 07/11/2024

### Introduction

This report documents the development of a SAFEM software that integrates MATLAB for core algorithm computation and FreeCAD for front-end modeling and visualization. The software is designed to streamline the process of finite element modeling, computation, and result visualization.

### Software Architecture

The software architecture consists of three main components:
1. **FreeCAD Front-End**: Utilized for model creation and obtaining essential FEA inputs. (Ongoing: specific modified UI for SAFEM, reference from "fcVM-workbench")
2. **MATLAB Core**: Handles the finite element algorithms and computations.
3. **Python Integration**: Facilitates the transfer of data between FreeCAD and MATLAB, as well as the visualization of results in FreeCAD. (Ongoing: the .py related to loading and Fourier transform function require the new UI inputs)

### Workflow

1. **Model Creation in FreeCAD**:
    - Models are created using FreeCAD, an open-source CAD software.
    - The model's dimensions, element division numbers, material parameters, and other necessary FEA inputs are defined.

2. **Data Transfer and Processing**:
    - A Python script is used to transfer the input data from FreeCAD to MATLAB.
    - The script processes the data to ensure compatibility with the MATLAB functions.

3. **Finite Element Computation in MATLAB**:
    - MATLAB scripts perform the core finite element computations.
    - The scripts handle various aspects of FEA such as meshing, material properties, boundary conditions, and solving the system equations.

4. **Result Visualization in FreeCAD**:
    - The computed results from MATLAB are transferred back to FreeCAD using the Python script.
    - FreeCAD is used to visualize the results, providing a comprehensive view of the analysis.

### MATLAB Scripts

The MATLAB scripts are organized based on their functionality, as listed below:

- **Dimension.m**: Defines the model dimensions.
- **Material.m**: Defines material properties.
- **Meshing.m**: Handles mesh generation.
- **Tire.m**: Specific to tire model analysis.
- **Temperature.m**: Manages temperature-related calculations.
- **Interface.m**: Defines interface properties.
- **Boundary.m**: Manages boundary conditions.
- **Contact_stress.m**: Computes the loading contact stress.
- **Fourier.m**: Performs Fourier transform analysis.
- **Process.m**: Handles the process parameters.
- **Solution.m**: Contains the main solution algorithm for FEA.
- **Results.m**: Manages the results data.
- **ResultsFigure.m**: Handles the generation of result figures.
- **test.m**: Used for initial testing and debugging.


### Conclusion

This integrated approach leveraging FreeCAD, MATLAB, and Python streamlines the finite element analysis process, from model creation to result visualization. The modularity and flexibility of the system allow for easy adaptation and expansion for various types of finite element problems. This development serves as a robust platform for conducting detailed and accurate finite element analyses.

### Further work
- **New UI**: design a modified FreeCAD UI for the special parameters input of SAFEM.
- **Python script**: achieve communication between FreeCAD and SAFEM. It retains all essential functions related to loading and result visualization, ensuring seamless interaction between these two software platforms.
