## FreeCAD GSoC 2024 Final development log: Implementation of Semi-Analytical Finite Element in FreeCAD for Pavement Infrastructure Modeling

## Kairen Shen 

## 08/08/2024
read the FEM workbench: https://wiki.freecad.org/FEM_Workbench
1. select the 2D FEM to present the pavement structure
2. use Contact Costraint to represent the interface between layers
3. material property use the elastic
issues:
1. extract element-node-coordinate information from FreeCAD
2. extract material, boundary, load information
 