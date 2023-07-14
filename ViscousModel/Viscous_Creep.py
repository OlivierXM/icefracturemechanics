"""
Dolfinx Implementation of PF Ice Fracture Mechanics

Using real load/displacement data, prescribe vertical traction at surface

Keys:
- Parallel support (HPC requires use of conda)
- Start of purely damage driven adaptive stepping (Fracture cases only)

strain = [XX, YY, XY]

"""
from mpi4py import MPI
from dolfinx import mesh, fem, nls, io, cpp, log, geometry
import numpy as np
import ufl
import petsc4py.PETSc as pets
import os
import sys
import datetime
import time
import math

# Custom Imports #
sys.path.append('../Utilities')
import Notes
import cdfx
import ViscousMats
import TestParameters

# MPI #
commWorld = MPI.COMM_WORLD
thisCore = commWorld.Get_rank()

notes = Notes.Notes()
fixes = Notes.Notes()
notes.addNotes("Modified Tolerance to 0.5e-5")
notes.addNotes("Specific Max Time of 4hr : Background Queue")

# Versioning #
appName = "Abaqus UMAT Implementation"
appVersion = "Version - 1.0 - DolfinX - Gc=1 Test Case"

# Saving #
refDir = ""
testDir = sys.argv[1]
subDir = sys.argv[2]
fullDir = refDir + subDir

# Mesh #
nCells = 25
# domain = mesh.create_rectangle(commWorld, [[0,0],[1,1]], [nCells,nCells], ghost_mode=mesh.GhostMode.shared_facet)
# domain2 = mesh.create_box(commWorld, [[1.25,0,0],[2.25,1,1]],[5,5,5], ghost_mode=mesh.GhostMode.shared_facet)
with io.XDMFFile(commWorld, "Bar.xdmf", "r") as xdmf:
    domain = xdmf.read_mesh(name="Grid")
    # cellTag = xdmf.read_meshtags(domain, name="Grid")

tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
f_to_c = domain.topology.connectivity(fdim, tdim)
num_cells = domain.topology.index_map(tdim).size_local + domain.topology.index_map(tdim).num_ghosts
num_facets = domain.topology.index_map(fdim).size_local + domain.topology.index_map(fdim).num_ghosts
totalCells = commWorld.allreduce(num_cells, op = MPI.SUM)
hmin = commWorld.allreduce(np.min(cpp.mesh.h(domain, tdim, range(num_cells))), op = MPI.MIN)
hmax = commWorld.allreduce(np.max(cpp.mesh.h(domain, tdim, range(num_cells))), op = MPI.MAX)
boxTree = cdfx.BoxTree(domain, commWorld) # Get x,y,z min/max coordinates of domain

# Files #
meshFile = io.XDMFFile(domain.comm, f"{fullDir}mesh.xdmf", "w")
meshFile.write_mesh(domain)
# with io.XDMFFile(domain.comm, f"{fullDir}meshref.xdmf", "w") as xdmf2:
#     xdmf2.write_mesh(domain2)
matFile = io.VTKFile(domain.comm, f"{fullDir}mats.pvd", "w")
facetFile = io.XDMFFile(domain.comm, f"{fullDir}facets.xdmf", "w")
cellFile = io.XDMFFile(domain.comm, f"{fullDir}cellTags.xdmf", "w")
damageFile = io.VTKFile(domain.comm, f"{fullDir}damageDir/damage.pvd", "w")
dispFile = io.VTKFile(domain.comm, f"{fullDir}dispDir/disp.pvd", "w")
sigmaFile = io.VTKFile(domain.comm, f"{fullDir}sigmaDir/sigma.pvd", "w")
elastFile = io.XDMFFile(domain.comm, f"{fullDir}elastDir/elastic.pvd", "w")
elastFile.write_mesh(domain)
hFile = io.VTKFile(domain.comm, f"{fullDir}hDir/h.pvd", "w")


logFileName = fullDir + "_LOGFILE.txt"
loadFileName = fullDir + "_LOADDISP.txt"
loadProfile = np.loadtxt("LoadFile.dat")
testParameters = TestParameters.TestParameters() # Custom test parameter class
fixes.addNotes("Finish vtk files")
fixes.addNotes("Anelastic response is weakly coupled to timestep")

# Files #
if (thisCore == 0):
     logFile = open(logFileName,'w')
     loadFile = open(loadFileName,'w')

# Constants #
in2m = 0.0254 # Convert inches to meters [m/in]
scalar1 = fem.Constant(domain, pets.ScalarType(1.0)) # Hold a fem.Constant of value 1

# Inputs #
l0 = testParameters._l0 # Nominal Crack Width [m] (25um)
gc_ice = testParameters._Gc # Critical Energy Release Rate [JDepth/m2]
gc_reinforced = (10**5)*gc_ice # Critical Energy Release Rate [J/m2]
gEff_ice = (2*l0)/gc_ice
gEff_reinforced = (2*l0)/gc_reinforced


# MaterialLibrary
globalTemp = 273-10 # -10C
matLib = ViscousMats.MaterialLibrary(globalTemp, domain, fromFile = testParameters)

# Viscous Parameters #
gamma = fem.Constant(domain, pets.ScalarType(matLib._gamma))
visc_constant = fem.Constant(domain, pets.ScalarType(matLib._viscous.ViscousConstant(globalTemp))) # 3.402639134892668e-13 if beta = 0.3

# Anelastic Parameters #
zeroMatrix = fem.Constant(domain, pets.ScalarType(( (0.0, 0.0),
                                                    (0.0, 0.0)
                                                    )))
zeroTens = fem.Constant(domain, pets.ScalarType(( (0.0, 0.0, 0.0),
                                                  (0.0, 0.0, 0.0),
                                                  (0.0, 0.0, 0.0),
                                                    )))
# Parameter Calculation

# Domain Generation #
Length = testParameters._LengthFactor * (boxTree._box[0, 1] - boxTree._box[0, 0]) # Length of bar (long dimension) [m]
Height = boxTree._box[1, 1] - boxTree._box[1, 0] # Height of bar (direction of 3-point loading) [m]
Depth = testParameters._Depth # width of bar (dimension into page) [m]
tol = 1e-12 # Tolerance for mesh tagging (not for computations) [m]

# Sim Conditions #
LoadingTol =(2 * 1.52E-05) # 3.5e-5
ReinforcedRadius = 0.002

ReinforcedRadiusSupport = 0.008
relSupportDistance = 1 * hmax
deltaLoad = 1.67e-5 # Displacement Rate [m/s]
dt1 = 1 # Starting timestep [s]
dt2 = dt1/4 # Secondary timestep [s]
dt3 = dt2/32 # Tertiary timestep [s] dt2/32
dt4 = dt3/1
staggerTime1 = 4e4 # 20
staggerTime2 = 7e4 # 25.75 /23
staggerTime3 = 9.0e4 # 37 / 25.75
totalTime = testParameters._TotalTime # Total Sim Time (larger than expected failure time) [s] 27.5
loadIteration_Step = 1 # Next Line to read
dt = dt1
dtPrevious = dt
# FunctionSpaces #
Fs_DG = fem.FunctionSpace(domain, ("DG", 0)) # Material Property Space
Fs_CG = fem.FunctionSpace(domain, ("CG", 1)) # Crack Solution Space (u_d)
Vs_CG = fem.VectorFunctionSpace(domain, ("CG", 1)) # Displacement Solution Space (u_k)
Vs_CG_3 = fem.VectorFunctionSpace(domain, ("DG", 0), dim = 3) # Displacement Solution Space (u_k)
Vs_Measure = fem.VectorFunctionSpace(domain, ("CG", 1), dim = 3)

# Functions #
u_d, u_k = fem.Function(Fs_CG), fem.Function(Vs_CG)
u_d_old, u_k_old = fem.Function(Fs_CG, name="Damage"), fem.Function(Vs_CG, name="Displacement")
v_d, v_k = ufl.TestFunction(Fs_CG), ufl.TestFunction(Vs_CG)
sigma_old, sigma_old_old, v_sigma = fem.Function(Vs_CG_3), fem.Function(Vs_CG_3), ufl.TestFunction(Vs_CG_3)
strain_old, v_strain = fem.Function(Vs_CG_3), ufl.TestFunction(Vs_CG_3)
sigma_delta, v_sigma_delta = fem.Function(Vs_CG_3), ufl.TestFunction(Vs_CG_3)
sigma_dev_delta, v_sigma_dev = fem.Function(Vs_CG_3), ufl.TestFunction(Vs_CG_3)
energy, energy_old, v_energy = fem.Function(Fs_DG, name="Energy"), fem.Function(Fs_DG), ufl.TestFunction(Fs_DG)
aEnergy = fem.Function(Fs_DG)
H_old, H_old_old, v_H_old = fem.Function(Fs_DG, name="H_old"), fem.Function(Fs_DG), ufl.TestFunction(Fs_DG)
u_k_old_sigma = fem.Function(Vs_CG)
u_d_old_sigma = fem.Function(Fs_CG)
anelasticStrainTotal = fem.Function(Vs_CG_3)
viscousStrainTotal, viscousStrainTotal_old = fem.Function(Vs_CG_3), fem.Function(Vs_CG_3)
elasticStrainTotal = fem.Function(Vs_CG_3)
stressMeasure = fem.Function(Vs_Measure)

# Pick this implementation or the other
matLib.Initialize(Vs_CG_3, Fs_CG, dt)

# Boundary Definitions #

class OmegaLoad(cdfx.SubDomain): # Reinforced load Areas
    def inside(self, x):
        x0 = x[0] - 0 # Loading is centered along length
        y0 = x[1] - Height # Loading is on top surface
        x1 = x[0] - 5*in2m # Support at Left End
        x2 = x[0] + 5*in2m # Support at Right End
        r1 = ReinforcedRadius
        r2 = ReinforcedRadiusSupport
        arg1 = (x0*x0+y0*y0 < r1*r1)
        arg2 = (x1*x1+x[1]*x[1] < r2*r2)
        arg3 = (x2*x2+x[1]*x[1] < r2*r2)
        arg4 = x[1] > 0.95*Height
        # return self.lor(self.lor(arg1,arg2),arg3)
        return self.lor(arg2,arg3)

class RigidConstraint(cdfx.SubDomain): # Bending Supports
    def inside(self, x):
        arg1 = self.near(x[1], 0, tol)
        arg2 = self.near(x[0],-0.5*Length + in2m, relSupportDistance)
        arg3  = self.near(x[0], 0.5*Length - in2m, relSupportDistance)
        return  self.land(arg1, self.lor(arg2, arg3))

class LoadingArea(cdfx.SubDomain): # Loading Area
    def inside(self, x):
        arg1 = self.near(x[1], Height, tol)
        arg2 = self.near(x[0], 0, LoadingTol)
        return  self.land(arg1, arg2)

class Crack(cdfx.SubDomain):
    def inside(self, x): # Initial Crack
        arg1 = self.near(x[0], 0, l0*2)
        arg2 = self.near(x[1], 0, l0*3)
        # arg1 = self.near(x[0], 0 - Length*0.05, l0*2)
        # arg2 = self.near(x[1], 0, 0.25*Height)
        return self.land(arg1, arg2)

class CenterLine(cdfx.SubDomain):
    def inside(self, x):
        return self.near(x[0], 0, tol)

class AllCheck(cdfx.SubDomain):
    def inside(self, x):
        arg1 = np.greater(abs(x[0]), -1)
        return arg1


notes.addNotes("Not currently enforcing fracture boundary conditions near crack")

# MatFunc()
# A helper function for defining gEff using lambda functions
# - x: The lambda delegate called inside gEff.interpolate()
def MatFunc(x):
    np_dist = np.sqrt(np.power(x[0] - 0, 2) + np.power(x[1] - Height, 2)) # np array of distances from the load point
    np_dist_log = np.less(np_dist, ReinforcedRadius).astype(np.int64) # a numpy int64 array of 1/0 determining if point is within ReinforcedRadius
    np_interior = gEff_reinforced # If the point is within the radius, assign it a constant value of gEff_reinforced
    np_exterior = gEff_ice + (gEff_reinforced-gEff_ice)*np.exp(-(1/ReinforcedRadius)*np.sqrt(np.abs(np.power(x[0] - 0, 2) + np.power(x[1] - Height, 2) - np.ones(x[0].shape)*ReinforcedRadius**2)))
    return np_dist_log*np_interior + (np.ones(x[0].shape) - np_dist_log)*np_exterior

# Materials #
gEff = cdfx.MaterialProperty(Fs_DG, default=gEff_ice)
gEff.interpolate(lambda x: MatFunc(x)) # Call the helper function to provide a smooth loading region
gEff.x.scatter_forward() # Share with MPI processes (lambda is only computed locally)
OmegaLoad().assign(gEff, gEff_reinforced) # override values at the supports
matFile.write([gEff._cpp_object])

# Cells #
cellFunction = cdfx.FacetFunction(domain, tdim, num_cells)
AllCheck().mark(cellFunction, 0)
cellTags = cellFunction.CreateMeshTag()
cellFile.write_mesh(domain)
cellFile.write_meshtags(cellTags)

# Facets #
## FacetFunction instance
# Hold the marked facets of topological dimension 2
facetFunction = cdfx.FacetFunction(domain, fdim, num_facets)
RigidConstraint().mark(facetFunction, 1)
CenterLine().mark(facetFunction, 4)
Crack().mark(facetFunction, 3)
LoadingArea().mark(facetFunction, 2)
facetTags = facetFunction.CreateMeshTag()
facetFile.write_mesh(domain)
facetFile.write_meshtags(facetTags)
facetFile.close()

# Setup dz and ds #

## ufl.Measure instance
# A facet integrator for integrating over tagged facets
dz = ufl.Measure('dx', domain = domain, subdomain_data = cellTags, metadata={"quadrature_degree": 5})
ds = ufl.Measure("ds", domain = domain, subdomain_data = facetTags)

# Load Expression #
loadArea = commWorld.allreduce(fem.assemble_scalar(fem.form(scalar1*ds(2))), op = MPI.SUM)
supportArea = commWorld.allreduce(fem.assemble_scalar(fem.form(scalar1*ds(1))), op = MPI.SUM)
totalVolume = commWorld.allreduce(fem.assemble_scalar(fem.form(scalar1*dz)), op = MPI.SUM)
displacementDriven = fem.Constant(domain, pets.ScalarType(0.0))
loadDriven = fem.Constant(domain, pets.ScalarType((0.0, 0.0)))
# Boundary Conditions #

## Displacement Boundary Condition List
# Hold the dirichlet boundary conditions for displacement
boundary_conditions_d = {1: {'Dirichlet': 1.0, 'Surface' : 3},
                        }

boundary_conditions_k = {1: {'Dirichlet': 0.0, 'Surface' : 1, 'Direction' : 1},
                         # 2: {'Dirichlet': displacementDriven, 'Surface' : 2, 'Direction' : 1},
                         3: {'Dirichlet': 0.0, 'Surface' : 2, 'Direction' : 0},
                        #  4: {'Dirichlet': 0.0, 'Surface' : 4, 'Direction' : 0},
                         }

## Displacement Boundary Condition Compiled List
# A python list used when declaring bcs in the solver
bcs_d = cdfx.DirichletBCs(boundary_conditions_d, Fs_CG, facetTags)
bcs_k = cdfx.DirichletBCs(boundary_conditions_k, Vs_CG, facetTags)

## Displacement Neumann Boundary Condition List
# A python list that when sum() can be added into the bilinear form
# integrals_N_k = [ufl.inner(loadDriven, v_k)*ds(2), -1*(loadArea/supportArea)*ufl.inner(loadDriven, v_k)*ds(1)]
integrals_N_k = [ufl.inner(loadDriven, v_k)*ds(2)]
integrals_N_d = []
# integrals_N_k = []
# Functions #

## epsilon()
# Return symmetric part of the gradient of u
# - u:a fem.Function from fem.VectorFunctionSpace
def epsilon(u):
    return 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

## Damage()
# (1-d)**2 with nonzero offset
# - d: a fem.function from fem.FunctionSpace
def Damage(d):
    return (1 - 1e-5)*(1.0 - d)**2 + 1e-5

## Damage_dd()
# Derivative of Damage()
# - d: a fem.function from fem.FunctionSpace
def Damage_dd(d):
    return (1 - 1e-5) * 2.0 * (d - 1.0)

## DevStress()
# Return the deviatoric part of the Cauchy stress
# - stressIn: the fem.Function from fem.TensorFunctionSpace
def DevStress(stressIn):
    vS = voigt2stress(stressIn)
    return vS - (1./3)*ufl.tr(vS)*ufl.Identity(tdim)

def DevStrain(u): # Deviatoric component of grad(u)
    localEps = epsilon(u)
    return localEps - (1./3) * ufl.tr(localEps)*ufl.Identity(tdim)


## Macaulay()
# Return the value if positive or negative
# - arg1: the value to evaluate
# - arg2: (1/-1), return value if +/- | else return 0
def Macaulay(arg1, arg2):
    return 0.5*(arg1 + arg2*abs(arg1))

## voigt3stress()
# Rewrite a vector into a symmetric tensor using Voigt notation
# - vec: the ufl object to rearrange (xx, yy, zz, xy, xz, yz)
def voigt3stress(vec):
    return ufl.as_tensor([[vec[0], vec[3], vec[4]],
                      [vec[3], vec[1], vec[5]],
                      [vec[4], vec[5], vec[2]]])

## voigt2stress()
# Rewrite a vector into a symmetric tensor using Voigt notation
# - vec: the ufl object to rearrange (xx, yy, xy)
def voigt2stress(vec):
    return ufl.as_tensor([[vec[0], vec[2]],
                      [vec[2], vec[1]]])

## ConvertToStress()
# Combine the effective strain and effect compliance tensor into effective stress using voigt3stress()
# - localStrain: The fem.Function from fem.TensorFunctionSpace containing a nxn strain matrix (xx, yy, xy)
# - complianceTensor: The return of ComplianceTensor(), a ufl.as_tensor
def ConvertToStress(localStrain, complianceTensor):
    localVec = ufl.as_vector([localStrain[0,0], localStrain[1,1], 2*localStrain[0,1]])
    return ufl.as_vector(ufl.dot(complianceTensor, localVec)) # dot product of (nxn) x (nx1) => nx1

## ViscousStrainChange()
# Return the stress dependent viscous strain change
# - stressIn: The stress at the previous state (fem.Function)
# - dtIn: The current time step [s]
def ViscousStrainChange(stressIn, dtIn):
    return dtIn * visc_constant * DevStress(stressIn)

## StrainHistory()
# Return the anelastic strain history ufl.as_tensor from ViscousMats.MaterialLibrary.GetStrain()
# - dtIn: The current time step [s]
def StrainHistory(dtIn):
    return matLib.GetStrain(dtIn, zeroMatrix)

## EffectiveStrain()
# Return the effective strain at the current timestep
# - u: The fem.Function to solve for
# - u_k_old: The fem.Function displacement from the previous time step
# - stressOld: The stress at the previous state (fem.Function)
# - dtIn: The current timestep [s]
def EffectiveStrain(u, u_k_old, stressOld, dtIn):
    return epsilon(u) - epsilon(u_k_old) - StrainHistory(dtIn) - ViscousStrainChange(stressOld, dtIn)

## sigTens()
# Return the damaged stiffness tensor when in compression
# When dIn = 0, unbroken, returns I * C_ijkl
# When dIn = 1, broken, returns volumetric component only
# - d: The fem.Function for the damage term
def sigCompr(dIn):
    d = Damage(dIn)
    _EV = matLib._E*matLib._anelastic._Vcoeff
    denom = (_EV*(d+2)-3*d+3)**2-(2*_EV*(d+2)+6*d+3)**2
    d11 = -(3*matLib._E*d*(d+2)*(2*_EV*(d+2)+6*d+3))/denom
    d12 = -(3*matLib._E*d*(d+2)*(_EV*(d+2)-3*d+3))/denom
    d33 = (matLib._E*d)/(_EV+1)
    damagedTensor = ufl.as_tensor([ [d11, d12, 0],
                                    [d12, d11, 0],
                                    [0,   0, d33],
                                    ])
    return damagedTensor * matLib._elastic._stiff_ufl

## sigTens()
# Return the damaged stiffness tensor when in tension or shear
# - d: The fem.Function for the damage term
def sigTens(d):
    return ufl.as_tensor(Damage(d) * matLib.GetCompliance())


## ComplianceTensor()
# Return the ufl.as_tensor by inverting the elastic+anelastic compliance tensor
def ComplianceTensor(u, d, _oldsigma):
    u_e = ufl.as_vector((1.0/matLib._E)*matLib._elastic._compl_ufl*_oldsigma) # Relying on old strain
    return ufl.conditional(ufl.lt(u_e[0]+u_e[1], 0.0), sigCompr(d), sigTens(d))

def ElasticComplianceTensor(u_k, u_d, _oldsigma):
    d = Damage(u_d)
    u_e = ufl.as_vector((1.0/matLib._E)*matLib._elastic._compl_ufl*_oldsigma) # Relying on old strain
    return ufl.conditional(ufl.lt(u_e[0]+u_e[1], 0.0), ElasticCompr(d), ElasticTens(d))

def ElasticTens(d):
    return ufl.as_tensor((1/(d * matLib._E)) * matLib._elastic._compl_ufl)

def ElasticCompr(d):
    denom = matLib._E * d * (d + 2)
    d11 = (2 * d + 1)/denom
    d12 = (d - 1)/denom
    d33 = 1 / (matLib._E * d)
    damagedTensor = ufl.as_tensor([ [d11, d12, 0],
                                    [d12, d11, 0],
                                    [0,   0, d33],
                                    ])
    return matLib._elastic._compl_ufl * damagedTensor

## Sigma()
# Return the updated stress tensor based on the previous stress state and the increase of stress
# - u: The fem.Function displacement to solve
# - u_k_old: The fem.Function displacement from the previous timestep
# - _oldsigma: The previous stress state tensor (fem.Function)
# - dtIn: The current timestep [s]
def Sigma(u_k, u_k_old, _oldsigma, u_d, u_d_old, dtIn):
    matLib.UpdateCompliance(dt)
    local_complianceNew = ComplianceTensor(u_k_old, u_d, _oldsigma) # Relying on old strain
    matLib.UpdateCompliance(dtPrevious)
    local_complianceOld = ufl.inv( ComplianceTensor(u_k_old, u_d_old, _oldsigma)) # on old strain
    matLib.UpdateCompliance(dt)
    return local_complianceNew * local_complianceOld * _oldsigma + ConvertToStress( EffectiveStrain(u_k, u_k_old, _oldsigma, dtIn), local_complianceNew)

def Sigma_old(u_k_old, _oldsigma, u_d, u_d_old, dtIn):
    matLib.UpdateCompliance(dt)
    local_complianceNew = ComplianceTensor(u_k_old, u_d, _oldsigma) # Relying on old strain
    matLib.UpdateCompliance(dtPrevious)
    local_complianceOld = ufl.inv( ComplianceTensor(u_k_old, u_d_old, _oldsigma)) # on old strain
    matLib.UpdateCompliance(dt)
    return local_complianceNew * local_complianceOld * _oldsigma

## SigmaVec()
# Return the (2n)x1 sigma vector as a ufl.as_vector
def SigmaVec(sigmaIn):
    return ufl.as_vector([sigmaIn[0,0],sigmaIn[1,1],sigmaIn[0,1]])

## Historic Energy Functional
# Enforce crack irreversibility by holding/updating total fracture energy
# - unew: The updated displacement as a fem.Function (vector)
# - Hold: The historic energy density (scalar)
# - energy: The viscous energy total (scalar)
def H(unew, Hold, energy):
    e_new = psi(energy)
    return ufl.conditional(ufl.lt(Hold, e_new), e_new, Hold)

def psi(vEnergy): # Fracture Potential
    elasticStrainLocal = elasticStrainTotal
    e1 = Macaulay(elasticStrainLocal[0] + elasticStrainLocal[1], 1)
    e2 = elasticStrainLocal - (1./3) * (elasticStrainLocal[0]+elasticStrainLocal[1]) * ufl.as_vector([1,1,0])
    posEnergy =  0.5 * matLib._Kv * (e1**2) + matLib._muinf * ufl.inner(e2, e2) + 0.0 * aEnergy
    return posEnergy + gamma*vEnergy

def ViscousEnergyChange(sigmaOld):
    avEnergy = 0.0
    return avEnergy + 0.5 * ufl.inner(SigmaVec(DevStress(sigmaOld)),(viscousStrainTotal - viscousStrainTotal_old))

def project(function, space):
    p = ufl.TrialFunction(space)
    q = ufl.TestFunction(space)
    a = ufl.inner(p, q) * dz
    L = ufl.inner(function, q) * dz
    problem = fem.petsc.LinearProblem(a, L, bcs = [])
    return problem.solve()

def dtFunc(d):
    return max((1.025/(1+math.exp((d-0.004)*1000)))**1.8, dt3)

# Problem Declaration #

# Displacement #
F_k = ufl.inner(voigt2stress(Sigma(u_k, u_k_old_sigma, sigma_old_old, u_d, u_d_old_sigma, dt)), epsilon(v_k))*dz - sum(integrals_N_k) #- ufl.inner(ufl.div(voigt2stress(Sigma_old(u_k_old_sigma, sigma_old_old, u_d, u_d_old_sigma, dt))),v_k)*dz

# Damage #
F_d = (u_d + gEff*Damage_dd(u_d)*H(u_k, H_old, energy))*v_d*dz + pets.ScalarType(4*l0*l0)*ufl.dot(ufl.grad(u_d),ufl.grad(v_d))*dz

# Historic Energy #
F_H_old = (H_old - H(u_k, H_old_old, energy))*v_H_old*dz

# Deviatoric #
F_dev = ufl.inner(sigma_dev_delta, v_sigma_dev)*dz - ufl.inner(SigmaVec(DevStress(sigma_delta)), v_sigma_dev)*dz

dtConst = fem.Constant(domain, pets.ScalarType(dtPrevious))
for i in range(matLib._M_size):
    j = 0
    localM = matLib._visc_solution_class[2*i+j]._m
    matLib._visc_solution_class[2*i+j]._vCoeff = fem.Constant(domain, pets.ScalarType(0.0))
    matLib._visc_solution_class[2*i+j]._vCoeff.value = matLib._anelastic.ReturnExpConst(localM, dtPrevious)[j]
    matLib._visc_solution_class[2*i + j].AssignF(
    ufl.inner(matLib._visc_solution_class[2*i + j]._u - matLib._visc_solution_class[2*i + j]._u_old *matLib._visc_solution_class[2*i+j]._vCoeff - \
    ((matLib._anelastic._Jd[localM]*matLib._anelastic._taud)/(dtConst*matLib._anelastic._Ld[localM]))*matLib._Dd*(1-matLib._visc_solution_class[2*i+j]._vCoeff)*ufl.dot(matLib._elastic._compl_ufl,sigma_dev_delta)
    ,matLib._visc_solution_class[2*i + j]._v) * dz, commWorld)
    j = 1
    localM = matLib._visc_solution_class[2*i+j]._m
    matLib._visc_solution_class[2*i+j]._vCoeff = fem.Constant(domain, pets.ScalarType(0.0))
    matLib._visc_solution_class[2*i+j]._vCoeff.value = matLib._anelastic.ReturnExpConst(localM, dtPrevious)[j]
    matLib._visc_solution_class[2*i + j].AssignF(
    ufl.dot(matLib._visc_solution_class[2*i + j]._u - matLib._visc_solution_class[2*i + j]._u_old *matLib._visc_solution_class[2*i+j]._vCoeff - \
    ((matLib._anelastic._Jgb[localM]*matLib._anelastic._taugb)/(dtConst*matLib._anelastic._Lgb[localM]))*matLib._Dgb*(1-matLib._visc_solution_class[2*i+j]._vCoeff)*ufl.dot(matLib._elastic._compl_ufl,sigma_dev_delta)
    ,matLib._visc_solution_class[2*i + j]._v) * dz, commWorld)


# Problem Creation #
Prob_k = fem.petsc.NonlinearProblem(F_k, u_k, bcs=bcs_k)
Prob_d = fem.petsc.NonlinearProblem(F_d, u_d, bcs=bcs_d)
Prob_dev = fem.petsc.NonlinearProblem(F_dev, sigma_dev_delta)
Prob_H_old = fem.petsc.NonlinearProblem(F_H_old, H_old)
# Solver Creation #

## Displacement #
Solver_k = nls.petsc.NewtonSolver(commWorld, Prob_k)
Solver_k.convergence_criterion = "residual"
Solver_k.atol = 1e-5
Solver_k.max_it = 25
Solver_k.report = True

# Damage #
Solver_d = nls.petsc.NewtonSolver(commWorld, Prob_d)
Solver_d.convergence_criterion = "residual"
Solver_d.atol = 1e-8
Solver_d.max_it = 25
Solver_d.report = True

# Deviatoric #
Solver_dev = nls.petsc.NewtonSolver(commWorld, Prob_dev)
Solver_dev.convergence_criterion = "incremental"
Solver_dev.atol = 1e-5
Solver_dev.max_it = 100
Solver_dev.report = False

# H_old #
Solver_H_old = nls.petsc.NewtonSolver(commWorld, Prob_H_old)
Solver_H_old.convergence_criterion = "incremental"
Solver_H_old.atol = 1e-8
Solver_H_old.report = False

# Write Log File #
if (thisCore == 0):
    logFile.write("{0}\n".format(appName))
    logFile.write("{0}\n".format(appVersion))
    logFile.write(str(datetime.datetime.now()) + "\n")
    logFile.write("Bar is {:.3f}in by {:.3f}in by {:.3f}in\n".format(Length/in2m,Height/in2m,Depth/in2m))
    logFile.write("Bar is {:.3f}mm by {:.3f}mm by {:.3f}mm\n".format(Length*1000,Height*1000,Depth*1000))
    logFile.write("Total number of processes is {0}".format(commWorld.Get_size()))
    logFile.write("\n-----Crack-----\n")
    logFile.write("\tCrack length, l0 = {:.2E}\n".format(l0))
    logFile.write("\tCritical energy rate, Gc_ice = {:.2f}\n".format(gc_ice))
    logFile.write("\tCritical energy rate, Gc_reinforced = {:.2f}\n".format(gc_reinforced))
    logFile.write("\n-----Constitutive-----\n")
    logFile.write("\tYoung's Modulus, E0 = {:.2E}\n".format(matLib._E))
    logFile.write("\tPoisson's ratio, nu0 = {:.2f}\n".format(matLib._nu))
    logFile.write("\tViscous contribution gamma, gamma = {:.2f}\n".format(matLib._gamma))
    logFile.write("\n-----Mesh-----\n")
    logFile.write("\tMin cell size = {:.2E}\n".format(hmin))
    logFile.write("\tTotal number of cells: {:.0f}\n".format(totalCells))
    logFile.write("\n-----Loading-----\n")
    logFile.write("\tNominal timestep, dt1 = {:.5f}\n".format(dt1))
    logFile.write("\tReduced timestep, dt2 = {:.5f}\n".format(dt2))
    logFile.write("\tReduced timestep, dt3 = {:.5f}\n".format(dt3))
    logFile.write("\tStagger Time1, t1 = {:.2f}\n".format(staggerTime1))
    logFile.write("\tStagger Time2, t2 = {:.2f}\n".format(staggerTime2))
    logFile.write("\tStagger Time3, t3 = {:.2f}\n".format(staggerTime3))
    logFile.write("\n-----Modifications-----\n")
    logFile.write("\tAnelastic Stiffness Mod, _JMod = {:.2f}\n".format(matLib._JMod))
    logFile.write("\tAnelastic Time Mod, _LMod = {:.2f}\n".format(matLib._LMod))
    logFile.write("\tViscous Mod, _BMod = {:.2f}\n".format(matLib._BMod))
    logFile.write("\n-----Notes-----\n")
    logFile.write(notes.logNotes())
    logFile.write("----------------\n")
    logFile.close()

    # Write Load File #
    loadFile.write("Time [s]\tDisplacement [m]\tActualLoad [N]\tStress [MPa]\tIterations [-]\tSim Time [s]\tStep Error\tTotal Damage\n")
    loadFile.close()

## Loop ##
t = 0
dt = dtFunc(0)
dtPrevious = dtFunc(0)
deltaLoad0 = 0
t0 = time.time()
itertol = 1e-5 # Staggered Tolerance
log.set_log_level(log.LogLevel.INFO)
damageValue = 0.0
loadIteration = 0

# Originate Crack #
Solver_d.solve(u_d)
fixes.addNotes("Only placing crack")

while ((t < totalTime) and (damageValue < 1.25e-02)):
    # Update Load
    # displacementDriven.value = deltaLoad0

    t = loadProfile[loadIteration,0]
    loadDriven.value = (0.0, (1.00/testParameters._LengthFactor)*loadProfile[loadIteration,1]/(Depth*loadArea))
    # loadDriven.value = (0.0, -20/(Depth*loadArea))
    displacementDriven.value = deltaLoad0
    loadIteration += loadIteration_Step
    dt = loadProfile[loadIteration,0] - t   

    dtConst.value = dtPrevious
    # Update Anelastic Strains
    aEnergy_func = 0.0
    for i in range(matLib._M_size):
        for j in range(2):
            matLib._visc_solution_class[2*i + j].Solve()
            localM = matLib._visc_solution_class[2*i+j]._m
            vCoeff = matLib._anelastic.GetVcoeff(dtPrevious, localM)[j]
            localStrain = ufl.as_vector((1 - matLib._anelastic.ReturnExpConst(localM, dtPrevious)[j])*matLib._visc_solution_class[2*i + j]._u_old)
            aEnergy_func = aEnergy_func + ufl.inner((1./vCoeff) * matLib._viscoStiff_ufl * localStrain, localStrain)

    if not(dt == dtPrevious): # Update time dependent coefficients if timestep changes
        matLib.UpdatePreviousVCoeff(dt)

    # Solve for damage + displacement
    iter = 0
    maxIter = 25
    err = 1



    if (matLib._M_size > 0):
        aEnergy = project(0.5*aEnergy_func, Fs_DG)

    # Update Viscous Strains
    viscousStrainTotal.x.array[:] = viscousStrainTotal_old.x.array + project(SigmaVec(ViscousStrainChange(sigma_old_old, dt)), Vs_CG_3).x.array
    viscousStrainTotal.x.scatter_forward()

    # Update Viscous Energy
    energy.x.array[:] = energy_old.x.array + project(ViscousEnergyChange(sigma_old_old), Fs_DG).x.array
    energy.x.scatter_forward()

    # <----- Stagger Here ----->
    while (err > itertol):
        iter += 1
        # Solve for damage
        # Solver_d.solve(u_d)
        # Solve for k
        Solver_k.solve(u_k)

        # Solve for the current sigma using the current and previous displacement
        sigma_old.interpolate(fem.Expression(Sigma(u_k, u_k_old_sigma, sigma_old_old, u_d, u_d_old_sigma, dt), Vs_CG_3.element.interpolation_points()))
        sigma_old.x.scatter_forward()

        # Update Elastic Strains
        elasticStrainTotal.x.array[:] = project(ufl.dot(ufl.as_tensor(matLib.GetElasticComplianceVol()), sigma_old), Vs_CG_3).x.array
        elasticStrainTotal.x.scatter_forward()
        anelasticStrainTotal.x.array[:] = project(SigmaVec(epsilon(u_k) - voigt2stress(viscousStrainTotal) - voigt2stress(elasticStrainTotal)), Vs_CG_3).x.array
        anelasticStrainTotal.x.scatter_forward()

        err_k = cdfx.errornorm(commWorld, (u_k - u_k_old)**2 * dz)
        err_d = cdfx.errornorm(commWorld, (u_d - u_d_old)**2 * dz)

        err = max(err_k, err_d)
        print(err_k, err_d)

        u_k_old.x.array[:] = u_k.x.array
        u_d_old.x.array[:] = u_d.x.array
        # Solver_H_old.solve(H_old) # Updates H_old
        errFinal = commWorld.allreduce(err, op = MPI.MAX) #MPI.MAX?

        if (errFinal < itertol and thisCore == 0):
            print(f'Iterations: {iter} Total time {t}')

        if (iter > maxIter):
            if (thisCore == 0):
                logFile = open(logFileName,'a')
                logFile.write("Maximum iterations ({:.0f}) reached at timestep {:.4f} with error {:1.3E}".format(maxIter,t, errFinal))
                logFile.close()

            raise Exception("Maximum iterations ({:.0f}) reached".format(maxIter), "Refine timestep at {:.4f} with error {:1.3E}".format(t, errFinal))


    sigma_delta.x.array[:] = sigma_old.x.array - sigma_old_old.x.array # Do not move after #Assignments
    Solver_dev.solve(sigma_dev_delta)
    # Assignments
    sigma_old_old.x.array[:] = sigma_old.x.array
    u_d_old_sigma.x.array[:] = u_d.x.array
    u_k_old_sigma.x.array[:] = u_k.x.array
    energy_old.x.array[:] = energy.x.array
    H_old_old.x.array[:] = H_old.x.array
    viscousStrainTotal_old.x.array[:] = viscousStrainTotal.x.array
    # Check if file should be saved
    if (round(t*1e2) % 1 == 0):
        dispFile.write([u_k_old._cpp_object], t)
        damageFile.write([u_d_old._cpp_object], t)
        sigmaFile.write([project(sigma_old_old[1],Fs_DG)._cpp_object], t)
        # stressMeasure.interpolate(fem.Expression(Sigma(u_k, u_k_old_sigma, sigma_old_old, u_d, u_d_old_sigma, dt), Vs_Measure.element.interpolation_points()))

        dispValue = commWorld.allreduce(fem.assemble_scalar(fem.form(u_k_old_sigma[1]*ds(2)))/loadArea, op = MPI.SUM)
        loadValue = testParameters._LengthFactor * commWorld.allreduce(Depth*fem.assemble_scalar(fem.form(sigma_old_old[1]*ds(2))), op = MPI.SUM)
        stressValue = testParameters._LengthFactor * commWorld.allreduce(fem.assemble_scalar(fem.form(sigma_old_old[1]*ds(2)))/loadArea, op = MPI.SUM)
        damageValue = commWorld.allreduce(fem.assemble_scalar(fem.form(u_d*dz))/totalVolume, op = MPI.SUM)
        viscTotal = commWorld.allreduce(fem.assemble_scalar(fem.form(viscousStrainTotal[1]*ds(2))), op = MPI.SUM)
        elasticTotal = commWorld.allreduce(fem.assemble_scalar(fem.form(elasticStrainTotal[1]*ds(2))), op = MPI.SUM)
        anelasticTotal = commWorld.allreduce(fem.assemble_scalar(fem.form(anelasticStrainTotal[1]*ds(2))), op = MPI.SUM)
        # strain1Value = commWorld.allreduce(fem.assemble_scalar(fem.form(matLib._visc_solution_class[0]._u_old[1]*ufl.dx)), op = MPI.SUM)
        # strain2Value = commWorld.allreduce(fem.assemble_scalar(fem.form(matLib._visc_solution_class[1]._u_old[1]*ufl.dx)), op = MPI.SUM)
        # viscousValue = commWorld.allreduce(fem.assemble_scalar(fem.form(energy1[1]*ufl.dx)), op = MPI.SUM)
        if (thisCore == 0):
            loadFile = open(loadFileName,'a')
            loadFile.write("{:.4f}".format(t)) # 1
            loadFile.write("\t{:1.3E}".format(dispValue))
            loadFile.write("\t{:1.3E}".format(loadValue))
            loadFile.write("\t{:1.3E}".format(stressValue)) # 4
            loadFile.write("\t{0}".format(iter)) # iter
            loadFile.write("\t{:.1f}".format(time.time()-t0))
            loadFile.write("\t{:1.3E}".format(errFinal)) # 7
            loadFile.write("\t{:1.3E}".format(damageValue))
            loadFile.write("\t{:1.3E}".format(viscTotal))
            loadFile.write("\t{:1.3E}".format(anelasticTotal))
            loadFile.write("\t{:1.3E}".format(elasticTotal))
            loadFile.write("\n")
            loadFile.close()
    cdfx.MPIprint(thisCore, 0, f"Time is {t}")
    dtPrevious = dt


## Editing ##
cdfx.MPIprint(0, thisCore, fixes.logNotes())
cdfx.MPIprint(0, thisCore, f"Total time is {commWorld.allreduce(time.time() - t0, op = MPI.MAX)}")
## End Script ##
