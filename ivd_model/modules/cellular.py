# modules/cellular.py

import numpy as np
from ..config.parameters import CellularParams as CP
from scipy.integrate import solve_ivp

class CellularModel:
    """
    Cell viability and matrix synthesis model
    Based on:
    1. Neidlinger-Wilke et al. (2014) - DOI: 10.1016/j.jbbm.2013.04.010
    2. Huang & Gu (2008) - DOI: 10.1016/j.jbiomech.2008.03.009
    """
    
    def __init__(self):
        self.params = CP
        self.reset_state()
        
    def reset_state(self):
        """Initialize cellular state variables"""
        self.state = {
            'cell_density': self.params.CELL_DENSITY,
            'cell_viability': 1.0,
            'collagen': 1.0,
            'proteoglycan': 1.0,
            'mmp_activity': 0.1,
            'timp_activity': 0.1
        }
        
    def calculate_cell_viability(self, oxygen, glucose, ph, mechanical_stress):
        """Calculate cell viability based on microenvironment"""
        # Critical thresholds
        O2_CRIT = 1.0  # kPa
        GLU_CRIT = 0.5  # mmol/L
        PH_MIN = 6.8
        STRESS_MAX = 0.3  # MPa
        
        # Calculate survival factors
        f_o2 = np.minimum(1.0, oxygen/O2_CRIT)
        f_glu = np.minimum(1.0, glucose/GLU_CRIT)
        f_ph = np.minimum(1.0, (ph - PH_MIN)/0.6)
        f_stress = np.maximum(0.0, 1.0 - mechanical_stress/STRESS_MAX)
        
        # Combined viability factor
        viability_factor = f_o2 * f_glu * f_ph * f_stress
        
        return viability_factor
    
    def update_matrix_synthesis(self, oxygen, mechanical_stress, dt):
        """Update matrix synthesis rates"""
        # Oxygen-dependent synthesis
        O2_opt = 5.0  # kPa, optimal oxygen tension
        f_o2 = np.exp(-(oxygen - O2_opt)**2/8.0)
        
        # Mechanical regulation
        stress_opt = 0.1  # MPa, optimal stress
        f_stress = np.exp(-(mechanical_stress - stress_opt)**2/0.02)
        
        # Matrix synthesis rates
        synthesis_factor = f_o2 * f_stress * self.state['cell_viability']
        
        d_collagen = (self.params.k_COL * synthesis_factor - 
                     self.params.k_DEG_COL * self.state['mmp_activity']) * dt
        
        d_proteoglycan = (self.params.k_PG * synthesis_factor - 
                         self.params.k_DEG_PG * self.state['mmp_activity']) * dt
        
        # Update state
        self.state['collagen'] += d_collagen
        self.state['proteoglycan'] += d_proteoglycan
        
        # Enzyme activities
        self.state['mmp_activity'] = (0.1 + 0.9 * mechanical_stress/0.3 * 
                                    (1 - self.state['timp_activity']))
        self.state['timp_activity'] = 0.1 * self.state['cell_viability']
        
        return self.state
    
    def update_state(self, oxygen, glucose, ph, mechanical_stress, dt):
        """Update all cellular state variables"""
        # Update cell viability
        viability = self.calculate_cell_viability(oxygen, glucose, ph, 
                                                mechanical_stress)
        self.state['cell_viability'] = viability
        
        # Update cell density
        death_rate = 0.1 * (1 - viability)
        self.state['cell_density'] *= (1 - death_rate * dt)
        
        # Update matrix synthesis
        self.update_matrix_synthesis(oxygen, mechanical_stress, dt)
        
        return self.state