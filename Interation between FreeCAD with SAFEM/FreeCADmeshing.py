import FreeCAD
import FreeCADGui
import Fem
import FemGui

# Open the FreeCAD document
doc = FreeCAD.openDocument("path/to/your/document.FCStd")

# Get the FEM mesh object
mesh_obj = doc.getObject("FEMMeshGmsh")  
mesh = mesh_obj.FemMesh

# Get node information
nodes = mesh.Nodes
print("Nodes:")
for node_id, node_coords in nodes.items():
    print(f"Node ID: {node_id}, Coordinates: {node_coords}")

# Get element information
elements = mesh.VolumeElements
print("\nElements:")
for elem_id, node_ids in elements.items():
    print(f"Element ID: {elem_id}, Node IDs: {node_ids}")

# If there are face elements
faces = mesh.FaceElements
print("\nFace Elements:")
for face_id, node_ids in faces.items():
    print(f"Face ID: {face_id}, Node IDs: {node_ids}")

# If there are edge elements
edges = mesh.EdgeElements
print("\nEdge Elements:")
for edge_id, node_ids in edges.items():
    print(f"Edge ID: {edge_id}, Node IDs: {node_ids}")

# Save mesh information to CSV files

import csv

# Save node information to a CSV file
with open("nodes.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Node ID", "X", "Y", "Z"])
    for node_id, node_coords in nodes.items():
        writer.writerow([node_id, *node_coords])

# Save element information to a CSV file
with open("elements.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Element ID", "Node IDs"])
    for elem_id, node_ids in elements.items():
        writer.writerow([elem_id, *node_ids])

# Save face element information to a CSV file
with open("face_elements.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Face ID", "Node IDs"])
    for face_id, node_ids in faces.items():
        writer.writerow([face_id, *node_ids])

# Save edge element information to a CSV file
with open("edge_elements.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Edge ID", "Node IDs"])
    for edge_id, node_ids in edges.items():
        writer.writerow([edge_id, *node_ids])
