import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
import meshio
from pathlib import Path
from ..config.parameters import MaterialParams

class SpatialModel:
    """
    Enhanced spatial model for IVD with distinct AF, NP, and CEP regions.
    Handles complex geometries from STL files and implements appropriate
    material properties and boundary conditions for each region.
    """
    
    def __init__(self):
        """Initialize the spatial model with real geometry"""
        self.regions = {}
        self.mesh = None
        self.material_params = MaterialParams()
        self.load_geometries()
        self.setup_fem()
        
    def load_geometries(self):
        """Load STL files for each region and merge them"""
        geometry_dir = Path(__file__).parent.parent / "geometry"
        
        # Load each region
        region_files = {
            'AF': geometry_dir / "af.stl",
            'NP': geometry_dir / "np.stl",
            'CEP': geometry_dir / "cep.stl"
        }
        
        meshes = {}
        for region, file_path in region_files.items():
            try:
                meshes[region] = meshio.read(str(file_path))
                print(f"Loaded {region} mesh: {len(meshes[region].points)} vertices, "
                      f"{len(meshes[region].cells[0].data)} elements")
            except Exception as e:
                print(f"Error loading {region} mesh: {str(e)}")
                continue
        
        # Merge meshes and keep track of element regions
        self.merge_meshes(meshes)
        
    def merge_meshes(self, meshes):
        """Merge individual region meshes into a single mesh while tracking regions"""
        all_points = []
        all_cells = []
        point_offset = 0
        self.element_regions = []  # Track which region each element belongs to
        
        for region, mesh in meshes.items():
            n_points = len(mesh.points)
            all_points.extend(mesh.points)
            
            # Adjust cell indices for the merged mesh
            region_cells = mesh.cells[0].data + point_offset
            all_cells.extend(region_cells)
            
            # Track element regions
            self.element_regions.extend([region] * len(mesh.cells[0].data))
            
            point_offset += n_points
        
        # Create merged mesh
        self.mesh = meshio.Mesh(
            points=np.array(all_points),
            cells=[("triangle", np.array(all_cells))]
        )
        
        # Store number of nodes and elements
        self.n_nodes = len(self.mesh.points)
        self.n_elements = len(all_cells)
        
        print(f"Merged mesh: {self.n_nodes} vertices, {self.n_elements} elements")
        
    def setup_fem(self):
        """Setup FEM matrices and initialize state variables"""
        self.nodes = self.mesh.points
        self.elements = self.mesh.cells[0].data
        
        # Initialize state variables
        self.state = {
            'displacement': np.zeros((self.n_nodes, 3)),  # 3D displacement field
            'stress': np.zeros((self.n_elements, 6)),     # σxx, σyy, σzz, σxy, σyz, σxz
            'strain': np.zeros((self.n_elements, 6)),     # εxx, εyy, εzz, γxy, γyz, γxz
            'pressure': np.zeros(self.n_nodes),
            'fluid_velocity': np.zeros((self.n_elements, 3))
        }
        
        # Setup region-specific material properties
        self.setup_material_properties()
        
        # Initialize boundary condition data
        self.identify_boundaries()
        
    def identify_boundaries(self):
        """Identify boundary nodes and surfaces"""
        # Get z coordinates
        z_coords = self.nodes[:, 2]
        max_z = np.max(z_coords)
        min_z = np.min(z_coords)
        tol = 0.1  # mm tolerance
        
        # Identify superior and inferior surfaces
        self.superior_nodes = np.where(np.abs(z_coords - max_z) < tol)[0]
        self.inferior_nodes = np.where(np.abs(z_coords - min_z) < tol)[0]
        
        # Identify outer AF surface nodes (for osmotic pressure boundary)
        self.identify_outer_af_surface()
        
    def identify_outer_af_surface(self):
        """Identify nodes on the outer surface of the AF"""
        # Find elements marked as AF
        af_elements = np.where(np.array(self.element_regions) == 'AF')[0]
        
        # Find boundary faces of AF elements
        boundary_faces = set()
        for elem in af_elements:
            faces = [
                tuple(sorted([self.elements[elem][i], self.elements[elem][(i+1)%3]]))
                for i in range(3)
            ]
            for face in faces:
                if face in boundary_faces:
                    boundary_faces.remove(face)
                else:
                    boundary_faces.add(face)
        
        # Get unique nodes from boundary faces
        self.outer_af_nodes = list(set([node for face in boundary_faces for node in face]))
        
    def apply_boundary_conditions(self, applied_load):
        """Apply boundary conditions for the mechanical problem"""
        ndof = 3 * self.n_nodes
        force = np.zeros(ndof)
        
        # Fixed bottom surface (inferior endplate)
        fixed_dofs = []
        for node in self.inferior_nodes:
            fixed_dofs.extend([3*node, 3*node+1, 3*node+2])
        
        # Applied load on superior surface
        for node in self.superior_nodes:
            force[3*node + 2] = -applied_load  # Negative for compression
        
        # Apply osmotic pressure on outer AF surface
        self.apply_osmotic_pressure(force)
        
        return np.array(fixed_dofs), force
        
    def apply_osmotic_pressure(self, force):
        """Apply osmotic pressure boundary conditions"""
        for node in self.outer_af_nodes:
            # Calculate normal vector at the node
            normal = self.calculate_normal_vector(node)
            
            # Calculate osmotic pressure magnitude
            pressure = self.calculate_osmotic_pressure(node)
            
            # Apply pressure force
            force[3*node:3*node+3] += pressure * normal
            
    def calculate_normal_vector(self, node):
        """Calculate normal vector at a surface node"""
        # Find connected elements
        connected_elements = [i for i, elem in enumerate(self.elements) 
                            if node in elem]
        
        # Calculate average normal from connected faces
        normal = np.zeros(3)
        for elem in connected_elements:
            # Get element vertices
            vertices = self.nodes[self.elements[elem]]
            
            # Calculate element normal
            v1 = vertices[1] - vertices[0]
            v2 = vertices[2] - vertices[0]
            elem_normal = np.cross(v1, v2)
            normal += elem_normal
            
        # Normalize
        return normal / np.linalg.norm(normal)
        
    def calculate_osmotic_pressure(self, node):
        """Calculate osmotic pressure based on fixed charge density"""
        # Get node position
        pos = self.nodes[node]
        
        # Base osmotic pressure (can be made more sophisticated)
        base_pressure = 0.15  # MPa
        
        # Vary with radial position
        r = np.sqrt(pos[0]**2 + pos[1]**2)
        r_max = np.max(np.sqrt(self.nodes[:,0]**2 + self.nodes[:,1]**2))
        
        return base_pressure * (1 - r/r_max)
        
    def update_stress_strain(self):
        """Update stress and strain state"""
        U = self.state['displacement'].flatten()
        
        for e in range(self.n_elements):
            # Get element nodes and displacements
            nodes = self.nodes[self.elements[e]]
            elem_dof = np.array([3*n + i for n in self.elements[e] for i in range(3)])
            u_e = U[elem_dof]
            
            # Calculate B matrix
            J = self.calculate_jacobian(nodes)
            dN = self.calculate_shape_derivatives(J)
            B = self.form_B_matrix(dN)
            
            # Calculate strain
            self.state['strain'][e] = B @ u_e
            
            # Calculate stress
            D = self.get_constitutive_matrix(e)
            self.state['stress'][e] = D @ self.state['strain'][e]
            
            # Update fluid velocity (Darcy's law)
            self.update_fluid_velocity(e)
            
    def get_constitutive_matrix(self, element):
        """Get constitutive matrix including fiber contribution if applicable"""
        E = self.material_props['E'][element]
        v = self.material_props['v'][element]
        
        # Get isotropic matrix
        D = self.isotropic_material_matrix(E, v)
        
        # Add fiber contribution for AF
        if self.element_regions[element] == 'AF':
            angle = self.material_props['fiber_angle'][element]
            density = self.material_props['fiber_density'][element]
            D_fiber = self.fiber_material_matrix(angle, density)
            D += D_fiber
            
        return D
        
    def update_fluid_velocity(self, element):
        """Update fluid velocity using Darcy's law"""
        # Get permeability
        k = self.material_props['k'][element]
        
        # Calculate pressure gradient
        nodes = self.elements[element]
        pressures = self.state['pressure'][nodes]
        
        # Get shape function derivatives
        J = self.calculate_jacobian(self.nodes[nodes])
        dN = self.calculate_shape_derivatives(J)
        
        # Calculate pressure gradient
        grad_p = dN.T @ pressures
        
        # Update fluid velocity (Darcy's law)
        self.state['fluid_velocity'][element] = -k * grad_p
        
    def solve_coupled_problem(self, applied_load, dt):
        """Solve coupled mechanical-fluid problem for one time step"""
        # Solve mechanical problem
        self.solve_mechanics(applied_load)
        
        # Update pressures based on volumetric strain
        self.update_pressures(dt)
        
    def update_pressures(self, dt):
        """Update pore pressures based on volumetric strain"""
        for e in range(self.n_elements):
            # Get volumetric strain
            εv = np.sum(self.state['strain'][e,:3])
            
            # Get material properties
            k = self.material_props['k'][e]
            β = self.material_props['β'][e]
            
            # Update pressure (simplified consolidation)
            nodes = self.elements[e]
            self.state['pressure'][nodes] += β * εv / (k * dt)
            
    def get_results(self):
        """Return current results for visualization or analysis"""
        return {
            'displacement': self.state['displacement'],
            'stress': self.state['stress'],
            'strain': self.state['strain'],
            'pressure': self.state['pressure'],
            'fluid_velocity': self.state['fluid_velocity'],
            'regions': self.element_regions
        }
