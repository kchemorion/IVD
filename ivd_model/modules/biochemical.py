# modules/biochemical.py

import numpy as np
from ..config.parameters import BiochemicalParams as BP
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve

class BiochemicalModel:
    def __init__(self, spatial_discretization=50):
        self.params = BP
        self.nx = spatial_discretization
        self.reset_state()
        
    def reset_state(self):
        """Initialize concentration fields"""
        self.state = {
            'oxygen': np.ones(self.nx) * 5.8,    # kPa
            'glucose': np.ones(self.nx) * 4.0,   # mmol/L
            'lactate': np.ones(self.nx) * 1.0,   # mmol/L
            'ph': np.ones(self.nx) * 7.1,         # pH units
            'water_content': np.ones(self.nx) * 0.8  # Initial water content (80%)
        }

    def create_diffusion_matrix(self, D, dx):
        """Create sparse matrix for diffusion equation with numerical stability"""
        main_diag = -2 * np.ones(self.nx)
        side_diag = np.ones(self.nx-1)
        
        # Create tridiagonal matrix
        A = diags([side_diag, main_diag, side_diag], [-1, 0, 1],
                shape=(self.nx, self.nx))
        
        # Safe calculation of diffusion coefficient
        dx = max(dx, 1e-6)  # Prevent too small dx
        diffusion_coeff = min(D/(dx*dx), 1e3)  # Prevent too large coefficients
        
        return csc_matrix(diffusion_coeff * A)

    def reaction_terms(self, c_o2, c_glu, c_lac, cell_density):
        """Calculate reaction terms for each species"""
        # Oxygen consumption (Michaelis-Menten kinetics)
        Km_o2 = 0.5  # Half-maximal concentration (kPa)
        R_o2 = -self.params.Q_O2 * cell_density * c_o2/(c_o2 + Km_o2)
        
        # Glucose consumption
        Km_glu = 0.3  # Half-maximal concentration (mmol/L)
        R_glu = -self.params.Q_GLU * cell_density * c_glu/(c_glu + Km_glu)
        
        # Lactate production (coupled to glucose consumption)
        R_lac = -2 * R_glu  # 2 lactate molecules per glucose
        
        return R_o2, R_glu, R_lac
        
    def update_concentrations(self, height, cell_density, dt):
        """Update concentration fields with water content regulation"""
        dx = max(height/self.nx, 1e-6)  # Prevent division by zero
        
        # Create diffusion matrices
        D_o2 = self.create_diffusion_matrix(self.params.D_O2, dx)
        D_glu = self.create_diffusion_matrix(self.params.D_GLU, dx)
        D_lac = self.create_diffusion_matrix(self.params.D_LAC, dx)
        
        # Calculate reaction terms
        R_o2, R_glu, R_lac = self.reaction_terms(
            self.state['oxygen'],
            self.state['glucose'],
            self.state['lactate'],
            cell_density
        )
        
        # Update species concentrations
        I = csc_matrix(np.eye(self.nx))
        
        # Safe update with bounds checking
        def safe_update(concentration, diffusion_matrix, reaction):
            new_conc = spsolve(I - dt*diffusion_matrix, concentration + dt*reaction)
            return np.clip(new_conc, 0, None)  # Prevent negative concentrations
        
        self.state['oxygen'] = safe_update(self.state['oxygen'], D_o2, R_o2)
        self.state['glucose'] = safe_update(self.state['glucose'], D_glu, R_glu)
        self.state['lactate'] = safe_update(self.state['lactate'], D_lac, R_lac)
        
        # Update pH based on lactate concentration
        self.state['ph'] = np.clip(7.4 - 0.1 * self.state['lactate'], 6.0, 7.4)
        
        # Update water content based on osmotic effects
        osmotic_pressure = -0.19 * self.state['lactate']  # Simplified osmotic pressure
        delta_water = -0.001 * osmotic_pressure * dt  # Water flux
        self.state['water_content'] = np.clip(
            self.state['water_content'] + delta_water,
            0.6,  # Minimum water content (60%)
            0.85  # Maximum water content (85%)
        )
        
        return self.state