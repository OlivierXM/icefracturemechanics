"""
    Simulation Test File
    Hold a set of constants that may be used by the Viscous Test Suite

"""
import numpy as np
## TestParameters()
# Copy any of the ViscousMats.MaterialConstants() that you want to modify here
# - Additionally add any other scalar constants
# - Some parameters are simulation specific (Required)  
class TestParameters():
    def __init__(self):
        # Optional #
        self._E = 6.50e9 # [Pa] 6.5
        self._nu = 0.33 # [-] 0.33
        self._JMod = 0.90 # [-] 0.90
        self._LMod = 0.81 # [-] 0.81
        self._BMod = 0.90 # [-] 0.90
        self._gamma = 0.0 # [-] 0.0
        # self._M = np.array([], dtype=np.int32) # Uncomment to remove Anelastic Mechanisms

        # Required #
        self._dt1 = 0.10**0.5556 # Starting timestep [s] 0.33**0.5556
        self._dt3 = 0.25/32 # Tertiary timestep [s] 0.25/32
        self._Gc = 4.50 # [J/m2] 4.5
        self._l0 = 50e-6 # [m] 50e-6
        self._Depth = 0.0254 # [m] 0.0254
        self._LoadRate = 1.67e-5 # [m/s] 1.67e-5
        self._LengthFactor = 1 # 1 if full bar, 2 if half bar
        self._TotalTime = 50 # [s] Simulation haults if t > TotalTime (avoid unnecessary sim time)
        self._primaryTol = 1e-5 # prior to dt3 timesteps (1e-5)
        self._secondaryTol = 1e-3 # post dt3 timesteps (1e-3)
        
