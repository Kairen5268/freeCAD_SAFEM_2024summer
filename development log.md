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
 
## 08/18/2024
1. read Introduction to Python: https://wiki.freecad.org/index.php?title=Introduction_to_Python
    There are several of ways to execute a Python program. But you can also execute it from the Python interpreter itself. For this, the interpreter must know where your program is. In FreeCAD the easiest way is to place your program in a folder that FreeCAD's Python interpreter knows by default, such as FreeCAD's user Mod folder:
    On Windows it is %APPDATA%\FreeCAD\Mod\, which is usually C:\Users\<username>\Appdata\Roaming\FreeCAD\Mod\.
    FreeCAD has a built-in Python interpreter. If you don't see the window labeled Python console as shown below, you can activate it under the View → Panels → Python console.

## 08/19/2024
1. read FreeCAD Scripting Basics: https://wiki.freecad.org/index.php?title=FreeCAD_Scripting_Basics
    For every App object in your document, there exists a corresponding Gui object. In fact the document itself has both an App and a Gui object. This, of course, only applies when you run FreeCAD with its full interface. In the command-line version no GUI exists, so only App objects are available. Note that the Gui part of objects is re-generated every time an App object is marked as 'to be recomputed' (for example when one of its parameters changes), so any changes made directly to the Gui object may be lost.

## 08/20/2024
Initially, search for the desired basic functionality through a script.
Analysis framework:
    a. develop the 3-D part and assembly the materials
    b. apply loading and boundary conditions
    c. conduct the meshing (by script) in the cross-section
    d. use above information for SAFEM analysis

1. read Part scripting: https://wiki.freecad.org/Part_scripting
    create the part by script
2. Code snippets: https://wiki.freecad.org/index.php?title=Code_snippets
    important!

## 08/21/2024
Develop the SAFEM by FreeCAD python script:

## 08/22/2024
Develop the SAFEM by FreeCAD python script:
    1. single layer structure

