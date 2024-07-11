import FreeCAD
import FreeCADGui
import Fem
import FemGui

# Open the FreeCAD document
doc = FreeCAD.openDocument("path/to/your/document.FCStd")

# Function to get material properties
def get_material_properties(obj):
    materials = obj.Material
    for mat in materials:
        print(f"Material Name: {mat.Name}")
        print(f"  Density: {mat['Density']}")
        print(f"  Young's Modulus: {mat['YoungsModulus']}")
        print(f"  Poisson's Ratio: {mat['PoissonRatio']}")

# Iterate through all objects in the document to find materials
for obj in doc.Objects:
    if obj.isDerivedFrom("Fem::FemMaterial"):
        print(f"Object Name: {obj.Name}")
        get_material_properties(obj)

# Save material properties to a CSV file
import csv

with open("material_properties.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Material Name", "Density", "Young's Modulus", "Poisson's Ratio"])
    for obj in doc.Objects:
        if obj.isDerivedFrom("Fem::FemMaterial"):
            materials = obj.Material
            for mat in materials:
                writer.writerow([mat.Name, mat['Density'], mat['YoungsModulus'], mat['PoissonRatio']])

print("Material properties have been extracted and saved to 'material_properties.csv'.")
