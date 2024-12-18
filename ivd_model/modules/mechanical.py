# modules/mechanical.py

import numpy as np
from scipy.integrate import solve_ivp
from ..config.parameters import MechanicalParams as MP

class MechanicalModel:
    """
    Biphasic mechanical model of the IVD
    Based on Malandrino et al. (2015)
    """
    def __init__(self):
        self.params = MP
        self.reset_state()
        
    def reset_state(self):
        """Initialize mechanical state variables"""
        self.state = {
            'height': self.params.DISC_HEIGHT,
            'strain': np.zeros(3),  # εxx, εyy, εzz
            'stress': np.zeros((3, 3)),  # Full stress tensor
            'pressure': 0.0,        # Fluid pressure
            'phi': self.params.phi0 # Current porosity
        }
        
    def calculate_stiffness_matrix(self):
        """Calculate the 3x3 stiffness matrix for normal strains"""
        E = self.params.E_AF  # Will be weighted average of AF/NP later
        v = self.params.v_AF
        
        # Form the 3x3 stiffness matrix for normal strains
        C = np.zeros((3, 3))
        factor = E/((1+v)*(1-2*v))
        
        # Fill the matrix
        for i in range(3):
            for j in range(3):
                if i == j:
                    C[i,j] = factor * (1-v)
                else:
                    C[i,j] = factor * v
                    
        return C
        
    def calculate_permeability(self):
        """
        Calculate strain-dependent permeability
        Based on exponential strain-dependence (Holmes-Mow model)
        """
        phi = self.state['phi']
        k0 = self.params.k0_AF  # Will be weighted average of AF/NP later
        M = 4.8  # Material parameter
        return k0 * np.exp(M * (phi - self.params.phi0))
        
    def calculate_osmotic_pressure(self):
        """
        Calculate osmotic pressure using van't Hoff equation
        """
        RT = 8.314 * 310  # Gas constant * Temperature
        c_fix = self.params.c_FCD * self.state['phi']
        return -0.19 * c_fix * RT
        
    def update_strain(self, applied_load, dt):
        """Update strain state with safety limits"""
        C = self.calculate_stiffness_matrix()
        k = self.calculate_permeability()
        p_osm = self.calculate_osmotic_pressure()
        
        def safe_deriv(t, y):
            strain = np.clip(y[:3], -0.3, 0.3)  # Limit strain range
            pressure = np.clip(y[3], -1e3, 1e3)  # Limit pressure range
            
            # Force vector (only axial load)
            force = np.array([0.0, 0.0, applied_load])
            force = np.clip(force, -1e3, 1e3)  # Limit force range
            
            # Strain rate with safety checks
            dstrain = np.zeros(3)
            try:
                dstrain = np.linalg.solve(C, force - pressure * np.ones(3))
                dstrain = np.clip(dstrain, -0.1, 0.1)  # Limit strain rate
            except np.linalg.LinAlgError:
                pass
            
            # Pressure rate with safety
            dp = np.clip(k * (applied_load - pressure - p_osm), -1e3, 1e3)
            
            return np.concatenate([dstrain, [dp]])
        
        # Initial conditions
        y0 = np.concatenate([
            np.clip(self.state['strain'], -0.3, 0.3),
            [np.clip(self.state['pressure'], -1e3, 1e3)]
        ])
        
        # Solve with tighter tolerances
        sol = solve_ivp(
            safe_deriv, [0, dt], y0,
            method='RK45',
            rtol=1e-6,
            atol=1e-8,
            max_step=dt/10
        )
        
        # Update state with bounds
        self.state['strain'] = np.clip(sol.y[:3,-1], -0.3, 0.3)
        self.state['pressure'] = np.clip(sol.y[3,-1], -1e3, 1e3)
        
        # Safe height update
        strain_z = np.clip(self.state['strain'][2], -0.3, 0.3)
        height_factor = np.clip(1 - strain_z, 0.5, 1.5)
        self.state['height'] = self.params.DISC_HEIGHT * height_factor
        
        # Update porosity with bounds
        vol_strain = np.sum(np.clip(self.state['strain'], -0.3, 0.3))
        denom = max(1 + vol_strain, 0.1)  # Prevent division by small numbers
        self.state['phi'] = np.clip(
            1 - (1 - self.params.phi0)/denom,
            0.3,  # Minimum porosity
            0.9   # Maximum porosity
        )
        
        return self.state
        
    def solve_mechanics(self, applied_load):
        """
        Solve the mechanical equilibrium problem for given load.
        
        Parameters:
        -----------
        applied_load : float
            Applied axial load in MPa
            
        Returns:
        --------
        dict
            Updated state dictionary containing:
            - height: current disc height (mm)
            - strain: strain components [εxx, εyy, εzz]
            - stress: stress components [σxx, σyy, σzz] (MPa)
            - pressure: fluid pressure (MPa)
            - phi: current porosity
        """        
        dt = 0.1  # Time step for mechanical equilibrium
        max_iter = 10
        tol = 1e-6
        
        for _ in range(max_iter):
            old_height = self.state['height']
            self.update_strain(applied_load, dt)
            
            # Update stresses - now store as 3x3 tensor
            C = self.calculate_stiffness_matrix()
            stress_elastic = np.dot(C, self.state['strain'])
            
            # Create full stress tensor
            self.state['stress'] = np.zeros((3, 3))
            for i in range(3):
                self.state['stress'][i,i] = stress_elastic[i] - self.state['pressure']
            
            # Add shear components if needed
            self.state['stress'][0,2] = self.state['stress'][2,0] = applied_load/2
            
            # Check convergence
            if np.abs(self.state['height'] - old_height) < tol:
                break
                
        return self.state