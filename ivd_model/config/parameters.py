class MaterialParams:
    """Material parameters for different IVD regions"""
    
    def __init__(self):
        # Annulus Fibrosus (AF)
        self.E_AF = 2.5          # Young's modulus (MPa)
        self.v_AF = 0.4          # Poisson's ratio
        self.k_AF = 1e-3         # Permeability (mm⁴/Ns)
        self.E_FIBER = 500.0     # Fiber modulus (MPa)
        
        # Nucleus Pulposus (NP)
        self.E_NP = 1.0          # Young's modulus (MPa)
        self.v_NP = 0.49         # Poisson's ratio (nearly incompressible)
        self.k_NP = 5e-4         # Permeability (mm⁴/Ns)
        
        # Cartilaginous Endplate (CEP)
        self.E_CEP = 5.0         # Young's modulus (MPa)
        self.v_CEP = 0.4         # Poisson's ratio
        self.k_CEP = 2e-3        # Permeability (mm⁴/Ns)
        
        # Fiber-related parameters
        self.FIBER_ANGLE_BASE = 30.0  # Base fiber angle in AF (degrees)
        self.FIBER_DISPERSION = 5.0   # Fiber angle dispersion (degrees)
        
        # Osmotic pressure parameters
        self.PI_0 = 0.15         # Base osmotic pressure (MPa)
        self.FCD_0 = 0.3         # Fixed charge density (mEq/mm³)
        
        # Biphasic parameters
        self.BIOT_COEF = 0.9     # Biot coefficient
        self.VOID_RATIO = 4.0    # Initial void ratio
        
        # Viscoelastic parameters
        self.TAU_1 = 3600        # Relaxation time constant 1 (s)
        self.TAU_2 = 60          # Relaxation time constant 2 (s)
        self.G_1 = 0.3           # Relative modulus 1
        self.G_2 = 0.1           # Relative modulus 2