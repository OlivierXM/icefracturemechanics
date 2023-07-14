"""
    Generate a dolfinx compatible xdmf/h5 from GMSH .geo file

    From https://github.com/OlivierXM/dolfinxTools
"""
import meshio
import os
import sys

## Offer default if no args passed
if (len(sys.argv) == 1):
    fileName = "Bar_Full_Fine"
    print(f"You have not provided an argument, the file {fileName}.geo will be meshed")
    userStop = input(f"Do you want to continue [y/n]:")
    if (userStop == "n" or userStop == "N"):
        print("Exiting: Aborted by User")
        quit()
else:
    fileName = sys.argv[1]

## Mesh with gmsh (2D)
dim = 2
os.system(f'gmsh {fileName}.geo -{dim} -format msh2')

## Read in mesh
msh = meshio.read(f'{fileName}.msh')

## Parse for 2D and 3D topology
for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data
    elif  cell.type == "tetra":
        tetra_cells = cell.data

## Scrub for physical markers
for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]
    elif key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]

## Write to file
if (dim == 2):
    # Strip msh.points of z with [:,:2]
    triangle_mesh = meshio.Mesh(points=msh.points[:,:2], cells=[("triangle", triangle_cells)])
    meshio.write("Bar.xdmf", triangle_mesh)
    print(f"Recommended HPC processes is : {int(round(len(triangle_cells)/40000,0))}!")

print(f"Output written to Bar.xdmf")