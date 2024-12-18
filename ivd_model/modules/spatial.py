# modules/spatial.py

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import meshio

class SpatialModel:
    """
    Spatial distribution model using FEM
    Based on:
    1. Malandrino et al. (2015) - DOI: 10.1016/j.jbiomech.2015.02.002
    2. Galbusera et al. (2011) - DOI: 10.1016/j.jbiomech.2011.04.021
    """
    
    def __init__(self, mesh_file=None, mesh_points=None):
        """
        Initialize spatial model
        
        Parameters:
        -----------
        mesh_file : str, optional
            Path to mesh file if using external mesh
        mesh_points : int, optional
            Number of points to use for default mesh (nx * ny)
        """
        if mesh_file:
            self.mesh = meshio.read(mesh_file)
        else:
            if mesh_points:
                # Calculate nx and ny to maintain aspect ratio
                ratio = 12/10  # height/radius ratio
                nx = int(np.sqrt(mesh_points/ratio))
                ny = int(nx * ratio)
            else:
                nx, ny = 20, 30
                
            self.create_default_mesh(nx, ny)
            
        self.setup_fem()
        
    def create_default_mesh(self, nx, ny):
        """Create a simple axisymmetric mesh"""
        x = np.linspace(0, 10, nx)  # radius
        y = np.linspace(0, 12, ny)  # height
        
        X, Y = np.meshgrid(x, y)
        points = np.column_stack((X.flatten(), Y.flatten()))
        
        # Create triangular elements
        cells = []
        for j in range(ny-1):
            for i in range(nx-1):
                n1 = j*nx + i
                n2 = j*nx + i + 1
                n3 = (j+1)*nx + i
                n4 = (j+1)*nx + i + 1
                cells.extend([[n1, n2, n3], [n2, n4, n3]])
                
        self.mesh = meshio.Mesh(
            points=points,
            cells={'triangle': np.array(cells)}
        )
        
    def setup_fem(self):
        """Setup FEM matrices"""
        self.nodes = self.mesh.points
        self.elements = self.mesh.cells_dict['triangle']
        
        self.n_nodes = len(self.nodes)
        self.n_elements = len(self.elements)
        
        # Initialize state variables at nodes
        self.state = {
            'displacement': np.zeros((self.n_nodes, 2)),
            'stress': np.zeros((self.n_elements, 3)),  # σrr, σzz, σrz
            'strain': np.zeros((self.n_elements, 3)),  # εrr, εzz, εrz
            'pressure': np.zeros(self.n_nodes),
            'oxygen': np.ones(self.n_nodes) * 5.8,
            'glucose': np.ones(self.n_nodes) * 4.0
        }
        
    def element_stiffness(self, nodes, E, v):
        """Calculate element stiffness matrix"""
        # Shape functions derivatives
        B = np.zeros((3, 6))
        
        # Calculate area and shape function derivatives
        x = nodes[:, 0]
        y = nodes[:, 1]
        
        area = 0.5 * abs(np.cross(nodes[1] - nodes[0], 
                                 nodes[2] - nodes[0]))
        
        # Form B matrix
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            B[0, 2*i] = (y[j] - y[k])/(2*area)
            B[1, 2*i+1] = (x[k] - x[j])/(2*area)
            B[2, 2*i] = B[1, 2*i+1]
            B[2, 2*i+1] = B[0, 2*i]
            
        # Material matrix
        D = np.zeros((3, 3))
        factor = E/((1+v)*(1-2*v))
        D[0,0] = D[1,1] = factor*(1-v)
        D[0,1] = D[1,0] = factor*v
        D[2,2] = factor*(1-2*v)/2
        
        return area * B.T @ D @ B
    
    def assemble_system(self, E, v):
        """Assemble global stiffness matrix"""
        ndof = 2 * self.n_nodes
        K = np.zeros((ndof, ndof))
        
        for e in range(self.n_elements):
            nodes = self.nodes[self.elements[e]]
            Ke = self.element_stiffness(nodes, E, v)
            
            # Assembly
            for i in range(3):
                for j in range(3):
                    rows = 2*self.elements[e][i] + np.arange(2)
                    cols = 2*self.elements[e][j] + np.arange(2)
                    K[np.ix_(rows, cols)] += Ke[2*i:2*i+2, 2*j:2*j+2]
                    
        return csr_matrix(K)
    
    def solve_mechanics(self, applied_load):
        """Solve mechanical problem"""
        # Material properties
        E = 2.5  # MPa
        v = 0.4
        
        # Assemble system
        K = self.assemble_system(E, v)
        
        # Apply boundary conditions
        # Bottom fixed
        fixed_nodes = np.where(self.nodes[:, 1] < 1e-6)[0]
        fixed_dofs = np.concatenate([2*fixed_nodes, 2*fixed_nodes + 1])
        
        # Top surface load
        top_nodes = np.where(self.nodes[:, 1] > np.max(self.nodes[:, 1]) - 1e-6)[0]
        force = np.zeros(2*self.n_nodes)
        force[2*top_nodes + 1] = -applied_load
        
        # Solve system
        free_dofs = np.setdiff1d(np.arange(2*self.n_nodes), fixed_dofs)
        U = np.zeros(2*self.n_nodes)
        U[free_dofs] = spsolve(K[np.ix_(free_dofs, free_dofs)], 
                              force[free_dofs])
        
        # Update state
        self.state['displacement'] = U.reshape(-1, 2)
        
        # Calculate strains and stresses
        self.update_stress_strain()
        
        return self.state
    
    def update_stress_strain(self):
        """Calculate element strains and stresses"""
        E = 2.5  # MPa
        v = 0.4
        
        for e in range(self.n_elements):
            nodes = self.nodes[self.elements[e]]
            U_e = self.state['displacement'][self.elements[e]].flatten()
            
            # Calculate B matrix
            B = np.zeros((3, 6))
            area = 0.5 * abs(np.cross(nodes[1] - nodes[0],nodes[2] - nodes[0]))
            
            for i in range(3):
                j = (i + 1) % 3
                k = (i + 2) % 3
                B[0, 2*i] = (nodes[j, 1] - nodes[k, 1])/(2*area)
                B[1, 2*i+1] = (nodes[k, 0] - nodes[j, 0])/(2*area)
                B[2, 2*i] = B[1, 2*i+1]
                B[2, 2*i+1] = B[0, 2*i]
            
            # Calculate strain
            self.state['strain'][e] = B @ U_e
            
            # Calculate stress
            D = np.zeros((3, 3))
            factor = E/((1+v)*(1-2*v))
            D[0,0] = D[1,1] = factor*(1-v)
            D[0,1] = D[1,0] = factor*v
            D[2,2] = factor*(1-2*v)/2
            
            self.state['stress'][e] = D @ self.state['strain'][e]