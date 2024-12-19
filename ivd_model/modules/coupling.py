import numpy as np
from scipy.interpolate import LinearNDInterpolator
from ..config.parameters import MaterialParams, BiochemicalParams

class ModelCoupling:
    """
    Handles coupling between different IVD model components:
    - Mechanical-Spatial: Stress/strain fields
    - Biochemical-Spatial: Nutrient transport
    - Cellular-Spatial: Mechanobiology
    """
    
    def __init__(self, spatial_model, mechanical_model, biochemical_model, cellular_model):
        self.spatial = spatial_model
        self.mechanical = mechanical_model
        self.biochemical = biochemical_model
        self.cellular = cellular_model
        
        self.material_params = MaterialParams()
        self.biochem_params = BiochemicalParams()
        
        # Initialize coupling matrices
        self.setup_coupling_matrices()
        
    def setup_coupling_matrices(self):
        """Setup interpolation matrices between different meshes"""
        # Get node positions from each model
        spatial_nodes = self.spatial.nodes
        mech_nodes = self.mechanical.nodes
        biochem_nodes = self.biochemical.nodes
        
        # Create KD-trees for efficient nearest neighbor search
        from scipy.spatial import KDTree
        self.spatial_tree = KDTree(spatial_nodes)
        self.mech_tree = KDTree(mech_nodes)
        self.biochem_tree = KDTree(biochem_nodes)
        
        # Store interpolation weights
        self.interp_weights = {}
        
        # Spatial to mechanical mapping
        spatial_to_mech = self.spatial_tree.query(mech_nodes, k=4)
        self.interp_weights['spatial_to_mech'] = {
            'indices': spatial_to_mech[1],
            'weights': self._calculate_weights(spatial_to_mech[0])
        }
        
        # Mechanical to spatial mapping
        mech_to_spatial = self.mech_tree.query(spatial_nodes, k=4)
        self.interp_weights['mech_to_spatial'] = {
            'indices': mech_to_spatial[1],
            'weights': self._calculate_weights(mech_to_spatial[0])
        }
        
        # Biochemical to spatial mapping
        biochem_to_spatial = self.biochem_tree.query(spatial_nodes, k=4)
        self.interp_weights['biochem_to_spatial'] = {
            'indices': biochem_to_spatial[1],
            'weights': self._calculate_weights(biochem_to_spatial[0])
        }
        
    def _calculate_weights(self, distances):
        """Calculate interpolation weights based on distances"""
        # Use inverse distance weighting
        weights = 1.0 / (distances + 1e-10)  # Avoid division by zero
        return weights / np.sum(weights, axis=1)[:, np.newaxis]
        
    def update_coupling(self, dt):
        """Update all coupled fields"""
        # 1. Update mechanical coupling
        self.update_mechanical_coupling()
        
        # 2. Update biochemical coupling
        self.update_biochemical_coupling(dt)
        
        # 3. Update cellular coupling
        self.update_cellular_coupling(dt)
        
    def update_mechanical_coupling(self):
        """Update mechanical-spatial coupling"""
        # Get mechanical state
        mech_stress = self.mechanical.state['stress']
        mech_strain = self.mechanical.state['strain']
        
        # Interpolate mechanical fields to spatial mesh
        spatial_stress = self._interpolate_field(
            mech_stress, 
            self.interp_weights['mech_to_spatial']
        )
        spatial_strain = self._interpolate_field(
            mech_strain,
            self.interp_weights['mech_to_spatial']
        )
        
        # Update spatial model state
        self.spatial.state['stress'] = spatial_stress
        self.spatial.state['strain'] = spatial_strain
        
        # Calculate fluid pressure from mechanical state
        self.update_fluid_pressure()
        
    def update_fluid_pressure(self):
        """Update fluid pressure based on mechanical state"""
        # Get volumetric strain
        strain = self.spatial.state['strain']
        vol_strain = np.sum(strain[:, :3], axis=1)
        
        # Get material properties
        k = np.array([self.material_params.k_AF if region == 'AF' else
                     self.material_params.k_NP if region == 'NP' else
                     self.material_params.k_CEP
                     for region in self.spatial.element_regions])
                     
        # Update pressure using Biot's theory
        alpha = self.material_params.BIOT_COEF  # Biot coefficient
        M = 2000.0  # Biot modulus (MPa)
        
        self.spatial.state['pressure'] = -alpha * M * vol_strain
        
    def update_biochemical_coupling(self, dt):
        """Update biochemical-spatial coupling"""
        # Get biochemical state
        concentrations = {
            'oxygen': self.biochemical.state['oxygen'],
            'glucose': self.biochemical.state['glucose'],
            'lactate': self.biochemical.state['lactate']
        }
        
        # Get fluid velocity from spatial model
        fluid_velocity = self.spatial.state['fluid_velocity']
        
        # Update advection-diffusion for each species
        for species in concentrations:
            # Get diffusion coefficient
            D = getattr(self.biochem_params, f'D_{species.upper()}_AF')
            
            # Solve advection-diffusion equation
            new_conc = self._solve_advection_diffusion(
                concentrations[species],
                fluid_velocity,
                D,
                dt
            )
            
            # Update biochemical model
            self.biochemical.state[species] = new_conc
            
    def _solve_advection_diffusion(self, concentration, velocity, diffusion, dt):
        """Solve advection-diffusion equation"""
        # Get mesh data
        nodes = self.spatial.nodes
        elements = self.spatial.elements
        
        # Calculate advection and diffusion terms
        n_nodes = len(nodes)
        dc_dt = np.zeros(n_nodes)
        
        for elem in range(len(elements)):
            # Get element nodes
            elem_nodes = elements[elem]
            
            # Calculate gradient of concentration
            grad_c = self._calculate_gradient(
                nodes[elem_nodes],
                concentration[elem_nodes]
            )
            
            # Advection term
            v = velocity[elem]
            adv_term = -np.dot(v, grad_c)
            
            # Diffusion term
            diff_term = diffusion * self._calculate_laplacian(
                nodes[elem_nodes],
                concentration[elem_nodes]
            )
            
            # Update rate of change
            dc_dt[elem_nodes] += (adv_term + diff_term) / len(elem_nodes)
        
        # Update concentration
        return concentration + dt * dc_dt
        
    def _calculate_gradient(self, nodes, values):
        """Calculate gradient using shape functions"""
        # Calculate Jacobian
        J = nodes[1:] - nodes[0]
        J_inv = np.linalg.inv(J)
        
        # Shape function derivatives
        dN = np.array([[-1, -1, -1],
                      [1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        
        # Transform to physical coordinates
        dN_dx = dN @ J_inv
        
        # Calculate gradient
        return dN_dx.T @ values
        
    def _calculate_laplacian(self, nodes, values):
        """Calculate Laplacian using shape functions"""
        # Simple approximation using finite differences
        return np.sum(values - values[0]) / (len(values) - 1)
        
    def update_cellular_coupling(self, dt):
        """Update cellular-spatial coupling"""
        # Get mechanical environment
        stress = self.spatial.state['stress']
        strain = self.spatial.state['strain']
        
        # Get biochemical environment
        oxygen = self.biochemical.state['oxygen']
        glucose = self.biochemical.state['glucose']
        ph = self.biochemical.state['ph']
        
        # Update cell density based on mechanical stimulus
        self.update_cell_density(stress, strain, dt)
        
        # Update cell metabolism based on local environment
        self.update_cell_metabolism(oxygen, glucose, ph)
        
        # Update ECM production/degradation
        self.update_ecm_remodeling(stress, strain, dt)
        
    def update_cell_density(self, stress, strain, dt):
        """Update cell density based on mechanical environment"""
        # Calculate mechanical stimulus
        von_mises = self._calculate_von_mises(stress)
        
        # Cell proliferation rate based on mechanical stimulus
        baseline_rate = 0.1  # per day
        mechanical_factor = np.clip(von_mises / 0.1, 0, 2)  # Normalized by 0.1 MPa
        
        # Update cell density
        density = self.cellular.state['cell_density']
        max_density = self.cellular.params.MAX_CELL_DENSITY
        
        growth_rate = baseline_rate * mechanical_factor * (1 - density/max_density)
        self.cellular.state['cell_density'] += growth_rate * dt
        
    def _calculate_von_mises(self, stress):
        """Calculate von Mises stress"""
        s11, s22, s33, s12, s23, s31 = stress.T
        return np.sqrt(0.5 * ((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2 +
                             6 * (s12**2 + s23**2 + s31**2)))
        
    def update_cell_metabolism(self, oxygen, glucose, ph):
        """Update cell metabolism based on local biochemical environment"""
        # Calculate metabolic rates based on local conditions
        o2_consumption = self._calculate_consumption_rate(
            oxygen,
            self.biochem_params.O2_CRITICAL,
            self.biochem_params.O2_CONSUMPTION_MAX
        )
        
        glucose_consumption = self._calculate_consumption_rate(
            glucose,
            self.biochem_params.GLUCOSE_CRITICAL,
            self.biochem_params.GLUCOSE_CONSUMPTION_MAX
        )
        
        # pH effect on metabolism
        ph_factor = np.clip((ph - 6.0) / (7.4 - 6.0), 0, 1)
        
        # Update cellular metabolic state
        self.cellular.state['metabolic_rates'] = {
            'oxygen': o2_consumption * ph_factor,
            'glucose': glucose_consumption * ph_factor
        }
        
    def _calculate_consumption_rate(self, concentration, critical, max_rate):
        """Calculate consumption rate using Michaelis-Menten kinetics"""
        Km = critical / 2  # Michaelis constant
        return max_rate * concentration / (Km + concentration)
        
    def update_ecm_remodeling(self, stress, strain, dt):
        """Update ECM remodeling based on mechanical state"""
        # Calculate mechanical stimulus
        von_mises = self._calculate_von_mises(stress)
        
        # Calculate ECM production/degradation rates
        production_rate = self._calculate_ecm_production(von_mises)
        degradation_rate = self._calculate_ecm_degradation(von_mises)
        
        # Update ECM content
        ecm = self.cellular.state['ecm_content']
        d_ecm = (production_rate - degradation_rate) * dt
        self.cellular.state['ecm_content'] = np.clip(ecm + d_ecm, 0, None)
        
    def _calculate_ecm_production(self, stimulus):
        """Calculate ECM production rate based on mechanical stimulus"""
        # Mechanically-induced production
        baseline = 0.1  # baseline production rate
        mechanical_factor = np.clip(stimulus / 0.05, 0, 2)  # normalized by 0.05 MPa
        return baseline * (1 + mechanical_factor)
        
    def _calculate_ecm_degradation(self, stimulus):
        """Calculate ECM degradation rate based on mechanical stimulus"""
        # Stress-induced degradation
        baseline = 0.05  # baseline degradation rate
        stress_factor = np.clip(stimulus / 0.2 - 1, 0, None)  # threshold at 0.2 MPa
        return baseline * (1 + stress_factor)
        
    def _interpolate_field(self, field, weights):
        """Interpolate field using precomputed weights"""
        indices = weights['indices']
        w = weights['weights']
        return np.sum(field[indices] * w[..., np.newaxis], axis=1)
