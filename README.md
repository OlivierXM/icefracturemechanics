# Ice Fracture Mechanics

Source code for the paper [An isotropic viscoelastic phase field fracture model for flexural loading
of freshwater columnar ice](https://doi.org/10.1016/j.commatsci.2023.112401).

## Description
A phase field fracture model for polycrystalline ice based on the constitutive model from [A viscoelastic integral formulation and numerical implementation of an isotropic constitutive model of saline ice](https://doi.org/10.1016/j.coldregions.2019.102983) by O'Connor et al. Model is implemented in the open source software [Dolfinx](https://fenicsproject.org/).

## Dolfinx Installation
Find details of installation at [Dolfinx's Github](https://github.com/FEniCS/dolfinx#installation)

## Usage
- Install dolfinx, gmsh and meshio
- Run `CreateMesh.sh` with corresponding *.geo file
- Run `Run_Viscous.py {1} {2} {3}`
    - Number of MPI processes
    - Directory to save data
    - Test Name
    - e.g. `Run_Viscous.py 6 . Viscous_Fracture.py`

## Support
Please direct questions to [@OlivierXM](https://github.com/OlivierXM). Many utilities maintained at [Olivier's Github](https://github.com/OlivierXM/dolfinxTools).

## License
This is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this code. If not, see http://www.gnu.org/licenses/.

We defer the reader to licensing associated with the following
- [GMSH](https://gmsh.info/#Documentation)
- [Dolfinx](https://fenicsproject.org/)
- [ParaView](https://www.paraview.org/resources/)
- [Doxygen](https://www.doxygen.nl/index.html)
- [SymPy](https://www.sympy.org/en/index.html)
- [Meshio](https://pypi.org/project/meshio/1.2.0/)
