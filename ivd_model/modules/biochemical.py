import numpy as np
from scipy.sparse import csr_matrix, linalg as sp_linalg
from scipy.sparse.linalg import spsolve
from ..config.parameters import BiochemicalParams as BP

class BiochemicalModel:
    """
    3D biochemical transport model for the IVD with reaction-diffusion equations
    """
    def __init__(self, geometry_model):
        self.params = BP
        self.geometry = geometry_model
        self.setup_transport()
        self.reset_state()
        
    def reset_state(self):
        """Initialize biochemical state variables"""
        self.state = {
            'oxygen': np.ones(self.n_nodes) * self.params.O2_INITIAL,
            'glucose': np.ones(self.n_nodes) * self.params.GLUCOSE_INITIAL,
            'lactate': np.zeros(self.n_nodes),
            'ph': np.ones(self.n_nodes) * self.params.PH_INITIAL,
            'cell_viability': np.ones(self.n_nodes)
        }
        
    def setup_transport(self):
        """Initialize transport matrices and boundary conditions"""
        # Get mesh data
        self.nodes = self.geometry.mesh.points
        self.elements = self.geometry.mesh.cells[0].data
        self.n_nodes = len(self.nodes)
        self.n_elements = len(self.elements)
        
        # Setup diffusion coefficients and reaction terms
        self.setup_coefficients()
        
        # Create FEM matrices
        self.create_fem_matrices()
        
        # Setup boundary conditions
        self.setup_boundary_conditions()
        
    def setup_coefficients(self):
        """Setup spatially varying coefficients"""
        # Initialize coefficient arrays
        self.D = {
            'oxygen': np.zeros(self.n_elements),
            'glucose': np.zeros(self.n_elements),
            'lactate': np.zeros(self.n_elements)
        }
        
        # Get AF/NP regions from geometry
        af_mask = self.geometry.get_af_mask()
        np_mask = ~af_mask
        
        # Set region-specific diffusion coefficients
        for species in ['oxygen', 'glucose', 'lactate']:
            D_af = getattr(self.params, f'D_{species.upper()}_AF')
            D_np = getattr(self.params, f'D_{species.upper()}_NP')
            self.D[species][af_mask] = D_af
            self.D[species][np_mask] = D_np
            
    def create_fem_matrices(self):
        """Create FEM matrices for transport equations"""
        # Create mass and stiffness matrices for each species
        self.M = {}  # Mass matrices
        self.K = {}  # Stiffness matrices
        
        for species in ['oxygen', 'glucose', 'lactate']:
            M, K = self.assemble_transport_matrices(self.D[species])
            self.M[species] = M
            self.K[species] = K
            
    def assemble_transport_matrices(self, diffusion_coef):
        """Assemble FEM matrices for transport equation"""
        # Initialize sparse matrix arrays
        n_dof = self.n_nodes
        nnz = self.n_elements * 16  # Approximate number of non-zeros
        
        rows = np.zeros(nnz, dtype=int)
        cols = np.zeros(nnz, dtype=int)
        mass_data = np.zeros(nnz)
        stiff_data = np.zeros(nnz)
        
        # Loop over elements
        idx = 0
        for el in range(self.n_elements):
            # Get element nodes
            nodes = self.elements[el]
            coords = self.nodes[nodes]
            
            # Calculate element matrices
            Me, Ke = self.element_matrices(coords, diffusion_coef[el])
            
            # Add to global matrices
            for i in range(4):
                for j in range(4):
                    rows[idx] = nodes[i]
                    cols[idx] = nodes[j]
                    mass_data[idx] = Me[i,j]
                    stiff_data[idx] = Ke[i,j]
                    idx += 1
                    
        # Create sparse matrices
        M = csr_matrix((mass_data[:idx], (rows[:idx], cols[:idx])), 
                      shape=(n_dof, n_dof))
        K = csr_matrix((stiff_data[:idx], (rows[:idx], cols[:idx])), 
                      shape=(n_dof, n_dof))
                      
        return M, K
        
    def element_matrices(self, coords, D):
        """Calculate element mass and stiffness matrices"""
        # Calculate Jacobian
        J = coords[1:] - coords[0]
        detJ = abs(np.linalg.det(J))
        Jinv = np.linalg.inv(J)
        
        # Shape function gradients
        dN = np.array([
            [-1, -1, -1],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])
        
        # Transform gradients to physical coordinates
        dNdx = dN @ Jinv
        
        # Mass matrix
        Me = detJ/20 * np.array([
            [2, 1, 1, 1],
            [1, 2, 1, 1],
            [1, 1, 2, 1],
            [1, 1, 1, 2]
        ])
        
        # Stiffness matrix
        Ke = D * detJ * (dNdx @ dNdx.T)
        
        return Me, Ke
        
    def setup_boundary_conditions(self):
        """Setup boundary conditions for transport"""
        # Find boundary nodes
        boundary = self.geometry.get_boundary_points()
        
        # Set Dirichlet BCs at blood supply
        self.bc_nodes = boundary['superior']
        self.bc_values = {
            'oxygen': self.params.O2_BLOOD,
            'glucose': self.params.GLUCOSE_BLOOD,
            'lactate': self.params.LACTATE_BLOOD
        }
        
    def solve_transport(self, dt, cell_consumption):
        """
        Solve coupled transport equations for one time step
        
        Parameters:
        -----------
        dt : float
            Time step size
        cell_consumption : dict
            Cell consumption/production rates for each species
        """
        # Update each species
        for species in ['oxygen', 'glucose', 'lactate']:
            # Get current concentration
            c = self.state[species]
            
            # Add consumption/production terms
            R = cell_consumption[species]
            
            # Create system matrix and RHS
            A = self.M[species]/dt + self.K[species]
            b = (self.M[species]/dt) @ c + R
            
            # Apply boundary conditions
            for node in self.bc_nodes:
                A[node, :] = 0
                A[node, node] = 1
                b[node] = self.bc_values[species]
                
            # Solve system
            c_new = spsolve(A, b)
            
            # Update state
            self.state[species] = c_new
            
        # Update pH based on lactate concentration
        self.update_ph()
        
        # Update cell viability based on metabolic state
        self.update_cell_viability()
        
        # Update geometry visualization
        self.update_geometry()
        
    def update_ph(self):
        """Update pH based on lactate concentration"""
        # Simple pH model based on lactate concentration
        # Could be made more sophisticated with proper acid-base chemistry
        base_ph = self.params.PH_INITIAL
        lactate_effect = 0.1 * self.state['lactate']
        self.state['ph'] = base_ph - lactate_effect
        np.clip(self.state['ph'], 6.0, 7.4, out=self.state['ph'])
        
    def update_cell_viability(self):
        """Update cell viability based on metabolic state"""
        # Define critical thresholds
        o2_crit = self.params.O2_CRITICAL
        glucose_crit = self.params.GLUCOSE_CRITICAL
        ph_crit_low = self.params.PH_CRITICAL_LOW
        ph_crit_high = self.params.PH_CRITICAL_HIGH
        
        # Calculate survival factors
        f_o2 = np.clip(self.state['oxygen'] / o2_crit, 0, 1)
        f_glucose = np.clip(self.state['glucose'] / glucose_crit, 0, 1)
        f_ph = np.clip((self.state['ph'] - ph_crit_low) / 
                      (ph_crit_high - ph_crit_low), 0, 1)
        
        # Update viability
        death_rate = (1 - f_o2 * f_glucose * f_ph) * self.params.DEATH_RATE
        self.state['cell_viability'] *= (1 - death_rate)
        np.clip(self.state['cell_viability'], 0, 1, out=self.state['cell_viability'])
        
    def update_geometry(self):
        """Update geometry model with current solution"""
        for field in self.state:
            self.geometry.interpolate_to_vertices(
                field,
                self.state[field],
                self.nodes
            )
