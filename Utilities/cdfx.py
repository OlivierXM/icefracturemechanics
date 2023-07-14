"""
Classical DolfinX (Legacy FEniCS)

Used to adjust user from verbose used in dolfin to DolfinX

For example the class SubDomain is recreated here to provide ease of use
similary to that of legacy FEniCS.

From https://github.com/OlivierXM/dolfinxTools
"""
from dolfinx import fem
from dolfinx import mesh
from mpi4py import MPI
import petsc4py.PETSc as pets
import numpy as np
import ufl



## MPIprint()
# Ensure only info from a specific core is printed
# - a: Typically the calling core [int]
# - b: The target core to respond [int]
# - c: The message [string]
def MPIprint(a, b, c):
    if (a == b):
        print(c)

## Unused
class ScalarExpression:
    def __init__(self, value):
        self.const = value

    def eval(self, x):
        # Added some spatial variation here. Expression is sin(t)*x
        return np.full(x.shape[1], self.const)

## Unused
class VectorExpression:
    def __init(self, value):
        self.const = value

    def eval(self, x):
        return

## SubDomain
# A dolfinx equivalent of dolfin.SubDomain
class SubDomain:
    ## SubDomain.__init__()
    # Default class constructor without args
    def __init__(self):
        return

    ## Unused, unable to implement at this time.
    def on_boundary(self):
        return False

    ## SubDomain.lor()
    # A shorthand for np.logical_or()
    # - a,b: logical np arrays
    def lor(self, a, b): # Numpy Logical Or
        return np.logical_or(a, b)

    ## SubDomain.land()
    # A shorthand for np.logical_and()
    # - a,b: logical np arrays
    def land(self, a, b): # Numpy logical And
        return np.logical_and(a, b)

    ## SubDomain.inside()
    # A placeholder function to be overridden by the class instance.
    # Use a combination of SubDomain.lor(), SubDomain.land(), and SubDomain.near() to determine marked facets
    # - x: The point coordinates, provided when called by other methods such as SubDomain.mark()
    def inside(self, x): # Place Holder
        return True

    ## Unused
    def insideBoundary(self, x):
        print(self.inside(x))
        print(self._boundary)
        print(f"Inside {self.inside(x).shape}")
        print(self._boundary.shape)
        return self.land(self.inside(x), self._boundary)

    ## SubDomain.near()
    # Same functionality as dolfin.near()
    # Return whether a given value is absolutely close to another point
    # - xIn: The point to check
    # - val: The reference point
    # - atol: The absolute tolerance to be considered "near"
    def near(self, xIn, val, atol): # Same functionality as dolfin.near()
        return np.isclose(xIn, val, rtol=0, atol=atol)

    def assign(self, matProp, value): # Return Numpy array of matching indices
        markedArray = mesh.locate_entities(matProp._domain, matProp._dim, self.inside)
        matProp.x.array[markedArray] = np.full(len(markedArray), value)
        return

    ## SubDomain.mark()
    # Given a facetFunction, mark appropriate facets with value as determined by SubDomain.inside()
    # - facetFunction: The cdfx.FacetFunction to mark
    # - value: The value to assign to all facets
    def mark(self, facetFunction, value):
        if (not(self.on_boundary())):
            markedArray = mesh.locate_entities(facetFunction._domain, facetFunction._fdim, self.inside)
        else:
            self._boundary = np.array(mesh.compute_boundary_facets(facetFunction._domain.topology))
            markedArray = mesh.locate_entities(facetFunction._domain, facetFunction._fdim, self.insideBoundary)

        facetFunction._indices.append(markedArray)
        facetFunction._markers[markedArray] = value
        return

## Alternative DirichletBC for marked domains
def DirichletBCs(argList, fs, facetTags):
    fdim = fs.mesh.topology.dim - 1
    bcs = []
    for i in argList:
        if 'Dirichlet' in argList[i]:
            if isinstance(argList[i]['Dirichlet'], fem.Constant):
                value = argList[i]['Dirichlet']
            else:
                value = fem.Constant(fs.mesh, pets.ScalarType(argList[i]['Dirichlet']))

            if 'Direction' in argList[i]:
                facets = np.array(facetTags.indices[facetTags.values==argList[i]['Surface']])
                left_dofs = fem.locate_dofs_topological(fs.sub(argList[i]['Direction']), fdim, facets)
                bcs.append(fem.dirichletbc(value, left_dofs, fs.sub(argList[i]['Direction'])))
            else:
                facets = np.array(facetTags.indices[facetTags.values==argList[i]['Surface']])
                left_dofs = fem.locate_dofs_topological(fs, fdim, facets)
                bcs.append(fem.dirichletbc(value, left_dofs, fs))

    return bcs

## Scalar Material Property (Renaming of fem.Function)
class MaterialProperty(fem.Function):
    ## MaterialProperty()
    # Construct a fem.Function with args
    # - fs: Provide the fem.FunctionSpace for this property
    # - default: Does the Material Property have an initial global value?
    def __init__(self, fs, default = None):
        super().__init__(fs)
        self._domain = fs.mesh
        self._dim = self._domain.topology.dim
        if not(default == None):
            self.assignAll(default)

    def assign(self):
        print("assign() is not directly supported by MaterialProperty, use assign from cdfx.SubDomain")        

    ## MaterialProperty.CellTagMark()
    # Given a MeshTags object defined over the domain, use tags to assign values
    # - region: The dolfinx.mesh.MestTags object
    # - marker: The corresponding tag to assign to
    # - value: The value to assign to the tagged region
    def CellTagMark(self, region, marker, value):
        markedArray = region.find(marker)
        self.x.array[markedArray] = np.full(len(markedArray), value)      
        
    ## MaterialProperty.assignAll()
    # Assign every cell in this collection the same value
    # - value: The value to assign
    def assignAll(self, value):
        self.x.array[:] = np.full(len(self.x.array), value)

## FacetFunction
# MeshFunction For Facets
class FacetFunction:
    ## FacetFunction.__init__()
    # Default constructor with args
    # - domain: The dolfinx.mesh entity
    # - fdim: The facet dimension for this function [int]
    def __init__(self, domain, fdim, numEntities):
        self._domain = domain
        self._fdim = fdim
        self._indices = []
        self._numEntities = numEntities
        # self._markers = []
        self._markers = np.zeros(self._numEntities, dtype=np.int32)
            
    def CellTagMark(self, region, marker, value):
        markedArray = region.find(marker)
        self._indices.append(markedArray)
        self._markers[markedArray] = value   

    ## FacetFunction.CreateMeshTag()
    # Call this after using SubDomain.mark() to generate a dolfinx.mesh.meshtags object
    def CreateMeshTag(self):
        return mesh.meshtags(self._domain, self._fdim, np.arange(self._numEntities, dtype=np.int32), self._markers)

## BoxTree
# Get min/max x,y,z coordinates of mesh
class BoxTree:
    ## BoxTree.__init__()
    # The default constructor with args
    # - domain: The dolfinx.mesh object to assess
    # - commWorld: The MPI.COMM_WORLD communicator object
    def __init__(self, domain, commWorld):
        self._domain = domain
        self._dim = domain.topology.dim
        self._box = np.ndarray(shape = (self._dim,2), dtype = float)
        for i in range(self._dim):
            self._box[i,0] = domain.comm.allreduce(np.min(domain.geometry.x[:, i]), op=MPI.MIN)
            self._box[i,1] = domain.comm.allreduce(np.max(domain.geometry.x[:, i]), op=MPI.MAX)

## errornorm()
# Calculate the L2 error
# - formArg: Provide the ufl form to assemble
def errornorm(comm, formArg):
    error_local = fem.assemble_scalar(fem.form(formArg))
    return np.sqrt(comm.allreduce(error_local, op=MPI.SUM))

## Function Variant
# Inherit dolfinx.Function for shorthand calls
# - fs: The fem.FunctionSpace for this Function
class FunctionX(fem.Function):
    def __init__(self, fs):
        super().__init__(fs)
        self._fs = fs

    def interp(self, expressIn, fs = None):
        if (fs == None):
            fs = self._fs
        self.interpolate(fem.Expression(expressIn, fs.element.interpolation_points()))

## End SCRIPT ##
