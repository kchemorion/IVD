import numpy as np
import meshio
from pathlib import Path

class GeometryModel:
    """Handles the 3D geometry of the IVD using real STL data"""
    
    def __init__(self):
        self.mesh = None
        self.centers = None
        self.vertex_values = {}  # Store field values at vertices
        self._load_geometry()
        
    def _load_geometry(self):
        """Load the STL geometry"""
        stl_path = Path(__file__).parent.parent / "geometry" / "ivd.stl"
        self.mesh = meshio.read(str(stl_path))
        
        # Calculate cell centers for field visualization
        self.centers = np.mean(
            self.mesh.points[self.mesh.cells[0].data],
            axis=1
        )
        
        # Initialize vertex values
        self._initialize_fields()
    
    def _initialize_fields(self):
        """Initialize field variables on the mesh"""
        # Create fields for various quantities
        fields = [
            'stress', 'strain', 'pressure', 'oxygen',
            'glucose', 'lactate', 'ph', 'cell_density'
        ]
        
        for field in fields:
            self.vertex_values[field] = np.zeros(len(self.mesh.points))
    
    def interpolate_to_vertices(self, field_name, values, positions):
        """
        Interpolate field values from computation points to mesh vertices
        using inverse distance weighting
        """
        from scipy.interpolate import LinearNDInterpolator
        
        # Create interpolator
        interp = LinearNDInterpolator(positions, values, fill_value=np.mean(values))
        
        # Interpolate to vertices
        self.vertex_values[field_name] = interp(self.mesh.points)
        
    def get_visualization_data(self):
        """
        Get mesh data in a format suitable for visualization
        """
        return {
            'vertices': self.mesh.points,
            'faces': self.mesh.cells[0].data,
            'vertex_values': self.vertex_values
        }
        
    def get_boundary_points(self):
        """
        Get points on the superior and inferior surfaces for boundary conditions
        """
        # Get z coordinates
        z_coords = self.mesh.points[:, 2]
        
        # Find superior and inferior surfaces with tolerance
        tol = 0.1  # mm
        max_z = np.max(z_coords)
        min_z = np.min(z_coords)
        
        superior_points = self.mesh.points[np.abs(z_coords - max_z) < tol]
        inferior_points = self.mesh.points[np.abs(z_coords - min_z) < tol]
        
        return {
            'superior': superior_points,
            'inferior': inferior_points
        }
