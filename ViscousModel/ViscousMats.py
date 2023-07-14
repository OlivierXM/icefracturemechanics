"""
Viscous Mats

Hold the constants used for the viscoelastic isotropic fracture of columnar ice (2D/3D)

Author: Olivier Montmayeur

Notes:
- For now just assume that temperature is isothermally considered
- Ignore crystallographic effects
- Does not support 1D problems...
- Currently 4 MaterialConstants() are generated (3 for each mechanism and 1 for MaterialLibrary())
-- They are identical in **kwargs `fromFile` but not linked (instances are not references of each other)

ToDo: 
- Test 3D usage
- Uncouple grain + dislocation pairings (create inidividual anelastic mechanisms at a time)

"""
import numpy as np
import math
from dolfinx import fem, nls
import petsc4py.PETSc as pets
import ufl
## MaterialConstants()
# A shared class that contains all the material inputs for this study
# - Only the Anelastic model has its own parameters, but they should not be changed
class MaterialConstants():
    ## MaterialConstants.__init__
    # Default constructor with optional args. Most if not all simulation parameters should be edited here.
    def __init__(self, fromFile = None):
        self._E = 4.0e9 # 4.0e9
        self._nu = 0.33 # 0.33
        self._JMod = 1.18 # 1.18
        self._LMod = 0.81 # 0.81
        self._BMod = 0.90 # 0.90
        self._gamma = 0.0 # 0.0
        self._M = np.array([0,1,2,3,4], dtype=np.int32) # Identify used constants
        self._dtOverride = False # Should the instantaneous compliance be time step independent (True/False)

        # Recursively iterate through any defined parameters above
        if not(fromFile == None):
            for a in vars(fromFile):
                for b in vars(self):
                    if (a==b):
                        vars(self)[a] = vars(fromFile)[a] # Copy Left

        # Constants ! Will override unless moved above !
        self._bd = (1/self._LMod) * 1.205e-9 # [Pa-s] 1.205e-9
        self._Qd = 0.55 # [eV] 0.55
        self._Qgb = 1.32 # [eV] 1.32
        self._bgb = (1/self._LMod) * 8.337e-28 # [eV] 8.337e-28
        self._Dgb = self._JMod * 3.0e-11 # [Pa^-1] 3.0e-11
        self._K = 0.07 # [Pa] 0.07
        self._rho = self._JMod * 1.29e9 # [m^-2] 1.29e9
        self._omega = (1.0 / math.pi) # [-] 1/pi
        self._boltzmann = 8.61733326e-5 # [eV/K]
        self._beta = self._BMod * (1/self._JMod) * (1/self._LMod) * 0.3 # [-] 0.3
        self._burg = 4.52e-10 # [m] 4.52e-10
        self._stressType = "Strain" # Stress/Strain for 2D Studies
        self.DerivedConstants()

    ## MaterialConstants.DerivedConstants()
    # Calculate material properties that derive from user inputs
    def DerivedConstants(self):
        self._Dd = (self._rho*self._omega*(self._burg**2))/self._K
        self._M_size = self._M.size
        self._Kv = self._E/(3*(1-2*self._nu))
        self._muinf = self._E/(2*(1+self._nu))


## MaterialLibrary()
# This class contains all the models used in this study
# - Only the anelastic model requires a specific contructor arg, which can be changed later
class MaterialLibrary(MaterialConstants):
    def __init__(self, tempIn, domain, **kwargs):
        super().__init__(**kwargs)
        self._tdim = domain.topology.dim
        if (self._tdim == 1):
           raise Exception("Viscous Mats Exception!", "Domains of topological dimension 1 are not supported!")

        self._elastic = ElasticMaterial(self._tdim, **kwargs)
        self._anelastic = AnelasticMaterial(tempIn, self._tdim, **kwargs)
        self._viscous = ViscousMaterial(**kwargs)
        self._viscoCompl = np.matmul(self._elastic._compl, self._anelastic._DDEVDSIG)
        self._viscoStiff_ufl = ufl.inv(ufl.as_tensor(self._viscoCompl))
        self._visc_solution_class = []
        self._elastic_const = fem.Constant(domain, pets.ScalarType(self._elastic._compl))
        self.DimensionalConstants()

    ## MaterialLibrary.DimensionalConstants()
    # Calculate derived constants from the instantiation of this class
    def DimensionalConstants(self):
        self._Ahat = CreateDDEVDSIG(ndi = self._tdim, ntens = 3*(self._tdim - 1))
        self._Ahat_m = np.eye(3*(self._tdim-1)) - self._Ahat

    ## MaterialLibrary.UpdateVcoeff_all()
    # Update the anelastic time dependent coefficients
    # - dtIn: The current time step [s]
    def UpdateVcoeff_all(self, dtIn):
        self._anelastic.UpdateVcoeff(dtIn)

    ## MaterialLibrary.UpdatePreviousVCoeff()
    # Update the previous time step time evolution
    # - dtIn: The previous time step [s]
    def UpdatePreviousVCoeff(self, dtIn):
        for i in range(len(self._visc_solution_class)):
            self._visc_solution_class[i]._vCoeff.value = self._anelastic.ReturnExpConst(self._visc_solution_class[i]._m, dtIn)[self._visc_solution_class[i]._index]

    ## MaterialLibrary.GetCompliance()
    # Return the current elastic+anelastic compliance tensor
    def GetCompliance(self):
        return self._totalCompl

    ## MaterialLibrary.GetElasticComplianceVol()
    # Return the elastic compliance tensor
    def GetElasticComplianceVol(self):
        return (1/self._E) * self._elastic._compl

    ## MaterialLibrary.GetAnelasticCompliance()
    # Return the anelastic compliance tensor
    def GetAnelasticCompliance(self):
        return self._anelastic._Vcoeff * self._viscoCompl

    ## MaterialLibrary.UpdateCompliance()
    # Update the time dependent anelastic compliance terms
    def UpdateCompliance(self, dtIn):
        self.UpdateVcoeff_all(dtIn) # Function of time only
        self._totalCompl = ufl.as_tensor(np.linalg.inv((1/self._E) * self._elastic._compl + self._anelastic._Vcoeff * self._viscoCompl))
        return

    ## MaterialLibrary.Initialize()
    # Initialize the anelastic components of the material
    # - Vs: The vector functionspace to use (fem.VectorFunctionSpace)
    # - fs: The scalar functionspace to use (fem.FunctionSpace)
    # - dtIn: The current time step [s]
    def Initialize(self, Vs, fs, dtIn):
        if not((self._M_size == 0)):
            for i in np.nditer(self._M): # Generate for each mechanism a disl. and grainb. function
                self._visc_solution_class.append(AnelasticFormulation(Vs, fs, i))
                self._visc_solution_class[-1]._index = 0
                self._visc_solution_class[-1]._timeFactor = ((self._anelastic._Jd[i]*self._anelastic._taud)/(dtIn*self._anelastic._Ld[i]))*(1.0-self._anelastic.ReturnExpConst(i, dtIn)[0])*self._Dd
                self._visc_solution_class.append(AnelasticFormulation(Vs, fs, i))
                self._visc_solution_class[-1]._index = 1
                self._visc_solution_class[-1]._timeFactor = ((self._anelastic._Jgb[i]*self._anelastic._taugb)/(dtIn*self._anelastic._Lgb[i]))*(1.0-self._anelastic.ReturnExpConst(i, dtIn)[1])*self._Dgb

        self.UpdateCompliance(dtIn)


    ## MaterialLibrary.GetStrain()
    # Return the anelastic strain history
    # - dtIn: The current time step [s]
    # - strainResult: a fem.Constant nxn zero matrix (placeholder)
    def GetStrain(self, dtIn, strainResult): # return nxn strain matrix
        if (self._M_size == 0):
            return strainResult # Return the nxn zero matrix if empty

        localResult = ufl.as_vector(np.zeros(shape=(int(3*(self._tdim-1)),)))
        j = 0
        for i in np.nditer(self._M):
            strainConst = self._anelastic.ReturnExpConst(i, dtIn) # d, gb
            localResult += (1-strainConst[0])*self._visc_solution_class[0+2*j]._u_old + (1-strainConst[1])*self._visc_solution_class[1+2*j]._u_old
            j += 1

        if (self._tdim == 3):
            return ufl.as_tensor([[localResult[0],0.5*localResult[3],0.5*localResult[4]],
                                [localResult[3],localResult[1],0.5*localResult[5]],
                                [0.5*localResult[4],0.5*localResult[5],localResult[2]],
                                ])
        else:
            return ufl.as_tensor([[localResult[0],0.5*localResult[2]],
                                [0.5*localResult[2],localResult[1]],
                                ])

## AnelasticMaterial()
# Contains the anelastic dislocation/grain boundary mechanisms
class AnelasticMaterial(MaterialConstants):
    def __init__(self, tempIn, tdim, **kwargs):
        super().__init__(**kwargs)
        self._Jd = np.ndarray((11,), dtype = np.double)
        self._Ld = np.ndarray((11,), dtype = np.double)
        self._Jgb = np.ndarray((11,), dtype = np.double)
        self._Lgb = np.ndarray((11,), dtype = np.double)
        self._taud = 0.0
        self._taudgb = 0.0
        self._Vcoeff = 0.0        
        self.Initialize(tdim)
        self.UpdateT(tempIn)

    ## AnelasticMaterial.Initialize()
    # Fill the anelastic mechanisms arrays with the fit values
    # - Additionally construct the deviatoric compliance matrix
    # - This should only be called once as these elements never change
    def Initialize(self, tdim):
        self._Jd[0] = 0.153051613823
        self._Ld[0] = 0.123946808159
        self._Jd[1] = 0.262856249885
        self._Ld[1] = 0.568624442765
        self._Jd[2] = 0.125295501586
        self._Ld[2] = 10.2032848417       
        self._Jd[3] = 0.04876601448
        self._Ld[3] = 68.0302680589
        self._Jd[4] = 0.275714252196
        self._Ld[4] = 2.26763749762
        self._Jgb[0] = 0.00504206655717
        self._Lgb[0] = 9.30861136337e-05  
        self._Jgb[1] = 0.0190291033313
        self._Lgb[1] = 0.00196781245785   
        self._Jgb[2] = 0.00938818210621
        self._Lgb[2] = 0.000164576765353  
        self._Jgb[3] = 0.0719182554952
        self._Lgb[3] = 0.0216204617582 
        self._Jgb[4] = 0.0289490376587
        self._Lgb[4] = 0.00265811200142 
          
        # self._Jd[0] = 0.153051613823
        # self._Ld[0] = 0.123946808159
        # self._Jd[1] = 0.262856249885
        # self._Ld[1] = 0.568624442765
        # self._Jd[2] = 0.125295501586
        # self._Ld[2] = 10.2032848417
        # self._Jd[3] = 0.0177870852505
        # self._Ld[3] = 710.263118465
        # self._Jd[4] = 0.04876601448
        # self._Ld[4] = 68.0302680589
        # self._Jd[5] = 0.275714252196
        # self._Ld[5] = 2.26763749762
        # self._Jd[6] = 0.00491210016165
        # self._Ld[6] = 22095.1526199
        # self._Jd[7] = 0.00938818210621
        # self._Ld[7] = 0.000164576765353
        # self._Jd[8] = 0.0719182554952
        # self._Ld[8] = 0.0216204617582
        # self._Jd[9] = 0.0289490376587
        # self._Ld[9] = 0.00265811200142
        # self._Jgb[0] = 0.0535360920955
        # self._Lgb[0] = 0.016851039015
        # self._Jgb[1] = 0.00504206655717
        # self._Lgb[1] = 9.30861136337e-05
        # self._Jgb[2] = 0.247381070943
        # self._Lgb[2] = 0.432597433957
        # self._Jgb[3] = 0.050441253476
        # self._Lgb[3] = 37.4939113511
        # self._Jgb[4] = 0.0173454650851
        # self._Lgb[4] = 339.328500832
        # self._Jgb[5] = 0.320488178311
        # self._Lgb[5] = 1.63411473906
        # self._Jgb[6] = 0.00441728078533
        # self._Lgb[6] = 8403.22064591
        # self._Jgb[7] = 0.155741127722
        # self._Lgb[7] = 6.12109095522
        # self._Jgb[8] = 0.126410130072
        # self._Lgb[8] = 0.0964994887392
        # self._Jgb[9] = 0.0190291033313
        # self._Lgb[9] = 0.00196781245785

        self._DDEVDSIG = CreateDDEVDSIG(ndi=tdim, ntens=int(3*(tdim-1)))

    ## AnelasticMaterial.UpdateT()
    # Updates the material temperature (currently isothermal)
    # - tempIn: The global temperature [K]
    def UpdateT(self, tempIn):
        self._taud = (self._bd/self._K)*math.exp(self._Qd/(self._boltzmann*tempIn))
        self._taugb = self._bgb*math.exp(self._Qgb/(self._boltzmann*tempIn))
        print(self._taud, self._taugb)
        return

    ## AnelasticMaterial.UpdateVcoeff()
    # Recalculate the material anelastic compliance matrix
    # - dtIn: The new timestep [s]
    def UpdateVcoeff(self, dtIn):
        # Assume UpdateT has been called #
        self._Vcoeff = 0.0
        if (self._M.size == 0): # Catch case where no anelastic components are active
            return

        if (self._dtOverride): # Prevent modification of timestep
            dtIn = 1.0

        for m in np.nditer(self._M): # Iterate through each of the selected anelastic mechanisms
            jm_d = self._Dd * self._Jd[m]
            lm_d = self._Ld[m] / self._taud
            jm_gb = self._Dgb * self._Jgb[m]
            lm_gb = self._Lgb[m] / self._taugb
            self._Vcoeff += jm_d*(1.0 - (1.0/(dtIn*lm_d))*(1 - math.exp(-dtIn*lm_d))) + \
                 jm_gb*(1.0 - 1.0/(dtIn*lm_gb)*(1.0 - math.exp(-lm_gb*dtIn)))
      
        return

    ## AnelasticMaterial.GetVcoeff()
    # Return a ufl vector with the anelastic compliance terms
    # dtIn: The current timestep [s]
    # m: The mechanism number to return
    def GetVcoeff(self, dtIn, m):
        if (self._dtOverride): # Prevent modification of timestep
            dtIn = 1.0

        jm_d = self._Dd * self._Jd[m]
        lm_d = self._Ld[m] / self._taud
        jm_gb = self._Dgb * self._Jgb[m]
        lm_gb = self._Lgb[m] / self._taugb
        return ufl.as_vector([jm_d*(1.0 - (1.0/(dtIn*lm_d))*(1 - math.exp(-dtIn*lm_d))), jm_gb*(1.0 - 1.0/(dtIn*lm_gb)*(1.0 - math.exp(-lm_gb*dtIn)))])

    ## AnelasticFormulation.ReturnExpConst()
    # Return the exp(dt*l/tau) for the dislocation/grainboundary mechanisms
    # - m: The mechanism # (int) to associate
    # - dtIn: The timestep [s]
    def ReturnExpConst(self, m, dtIn):
        return np.array([(math.exp((-dtIn*self._Ld[m])/self._taud)),(math.exp((-dtIn*self._Lgb[m])/self._taugb))], dtype=np.float64)

## ViscousMaterial()
# Holds the viscous evolution constants based on MaterialConstants()
class ViscousMaterial(MaterialConstants):
    ## ViscousMaterial.__init__()self
    # Default argless constructor
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    ## ViscousMaterial.ViscousConstant()
    # Return the updated viscous constant based on the current temperature
    # - tempIn: The new temperature [K]
    def ViscousConstant(self, tempIn):
        return ((self._beta*self._rho*(self._omega**1.5)*(self._burg**2))/(self._bd))*math.exp(-(self._Qd/(self._boltzmann*tempIn)))

## ElasticMaterial()
# Holds the elastic compliance tensor
class ElasticMaterial(MaterialConstants):
    ## ElasticMaterial.__init__()
    # Default contstructor argless that generates the elastic compliance tensor
    # - This should be updated once temperature varying properties are included
    def __init__(self, tdim, **kwargs):
        super().__init__(**kwargs)
        if (tdim == 3):
            self._compl = np.array([[1, -self._nu, -self._nu, 0., 0., 0.],
                           [-self._nu, 1, -self._nu, 0., 0., 0.],
                           [-self._nu, -self._nu, 1, 0., 0., 0.],
                           [0., 0., 0., 2*(1+self._nu), 0., 0],
                           [0., 0., 0., 0., 2*(1+self._nu), 0.],
                           [0., 0., 0., 0., 0., 2*(1+self._nu)]])
        else:
            if (self._stressType == "Stress"):
                self._compl = np.array([[1, -self._nu, 0.],
                                       [-self._nu, 1, 0.],
                                       [0., 0., 2*(1+self._nu)]])
            else:
                self._compl = np.array([[1-self._nu**2, -self._nu*(self._nu+1), 0.],
                                       [-self._nu*(self._nu+1), 1-self._nu**2, 0.],
                                       [0., 0., 2*(1+self._nu)]])


        self._compl_ufl = ufl.as_tensor(self._compl)
        self._stiff_ufl = ufl.inv(self._compl_ufl)

## AnelasticFormulation()
# A custom class to handle each coupled dislocation/grainboundary anelastic mechanisms
# This class solves for the time evolution of the anelastic strain histories
# Each mechanism, potentially up to 20, will have its own:
# - fem.Function(), ufl.TestFunction()
# - fem.petsc.NonlinearProblem()
# - nls.petsc.NewtonSolver()
class AnelasticFormulation():
    ## AnelasticFormulation.__init__()
    # A constructor with args
    # - vs: The fem.VectorFunctionSpace to use. (n-1)*3 x 1
    # - fs: The fem.FunctionSpace to use.
    # - m: The mechanism # (int) to associate
    def __init__(self, vs, fs, m):
        self._u = fem.Function(vs)
        self._u_old = fem.Function(vs)
        self._v = ufl.TestFunction(vs)
        self._viscstrain = fem.Function(vs)
        self._viscstrainChange = fem.Function(vs)
        self._m = m
        self._timeFactor = 0.0
        self._index = 0
        self._problem = None
        self._solver = None
        self._vCoeff = 0.0
        self._vs = vs

    ## AnelasticFormulation.AssignF()
    # Assign the F=0 bilinear form for fem.petsc.NonlinearProblem()
    # - f: The blinear form
    # - commWorld: The MPI communicator
    def AssignF(self, f, commWorld):
        self._F = f
        self._problem = fem.petsc.NonlinearProblem(self._F, self._u)
        self._solver = nls.petsc.NewtonSolver(commWorld, self._problem)
        self._solver.convergence_criterion = "incremental"
        self._solver.atol = 1e-5
        self._solver.max_it = 100
        self._solver.report = False

    ## AnelasticFormulation.Solve()
    # Update the anelastic history
    def Solve(self):
        self._solver.solve(self._u)
        self._viscstrainChange.interpolate(fem.Expression((1.0 - self._vCoeff) * self._u_old, self._vs.element.interpolation_points()))
        # self._viscstrain.x.array[:] = self._viscstrain.x.array + #self._viscstrainChange.x.array
        # self._viscstrain.x.scatter_forward()
        self._u_old.x.array[:] = self._u.x.array


## CreateDDEVDSIG()
# Generate the deviatoric compliance matrix
# - ndi: for 2D ndi = 2, for 3D ndi = 3 (problem dimension )
# - ntens: for 2D ndi = 3, for 3D ndi = 6 (size of symmetric compliance tensor, [xx, yy, xy] for 2D ....)
def CreateDDEVDSIG(ndi=3, ntens=6):
    NDI = ndi
    NTENS = ntens
    a = np.ndarray((NTENS,NTENS), dtype=np.double)
    a[:,:] = 0.0
    for i in range(NDI):
        for j in range(NDI):
            a[i,j] = -1.0/3

        a[i,i] = 2.0/3

    for i in range(NDI, NTENS):
        a[i,i] = 1.0

    return a
