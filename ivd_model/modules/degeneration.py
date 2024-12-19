import numpy as np
from scipy.integrate import solve_ivp
from ..config.parameters import DegenerationParams as DP

class DegenerationModel:
    """
    Comprehensive degeneration model for IVD including:
    - Age-related changes
    - Injury response
    - Metabolic stress effects
    """
    def __init__(self, geometry_model):
        self.params = DP
        self.geometry = geometry_model
        self.setup_degeneration_model()
        
    def setup_degeneration_model(self):
        """Initialize degeneration model parameters"""
        # Get mesh information
        self.nodes = self.geometry.mesh.points
        self.elements = self.geometry.mesh.cells[0].data
        self.n_nodes = len(self.nodes)
        
        # Initialize degeneration state variables
        self.state = {
            # Tissue composition
            'gag_content': np.ones(self.n_nodes),  # Normalized GAG content
            'collagen_integrity': np.ones(self.n_nodes),  # Collagen network integrity
            'water_content': np.ones(self.n_nodes) * 0.8,  # Initial water content
            
            # Cell population
            'senescent_cells': np.zeros(self.n_nodes),  # Fraction of senescent cells
            'stem_cells': np.zeros(self.n_nodes),  # Resident stem cells
            
            # Inflammatory factors
            'mmps': np.zeros(self.n_nodes),  # Matrix metalloproteinases
            'cytokines': np.zeros(self.n_nodes),  # Inflammatory cytokines
            'growth_factors': np.zeros(self.n_nodes),  # Anabolic growth factors
            
            # Damage accumulation
            'mechanical_damage': np.zeros(self.n_nodes),  # Accumulated mechanical damage
            'oxidative_stress': np.zeros(self.n_nodes),  # Oxidative stress level
            
            # Degeneration grade
            'pfirrmann_grade': np.ones(self.n_nodes) * 1  # Initial healthy state
        }
        
    def simulate_age_related_degeneration(self, current_age, time_span):
        """
        Simulate age-related degeneration
        
        Parameters:
        -----------
        current_age : float
            Current age in years
        time_span : float
            Duration of simulation in years
        """
        # Age-related parameter changes
        def age_effects(age):
            return {
                'cell_death_rate': self.params.BASE_CELL_DEATH_RATE * 
                                 (1 + self.params.AGE_DEATH_FACTOR * (age - 20)),
                'gag_synthesis': self.params.BASE_GAG_SYNTHESIS * 
                               np.exp(-self.params.AGE_GAG_DECAY * (age - 20)),
                'senescence_rate': self.params.BASE_SENESCENCE_RATE * 
                                 (1 + self.params.AGE_SENESCENCE_FACTOR * (age - 20))
            }
            
        # Time points for simulation
        times = np.linspace(0, time_span * 365, int(time_span * 12))  # Monthly points
        
        for t in times:
            current_time = current_age + t/365
            effects = age_effects(current_time)
            
            # Update cellular state
            self.update_cellular_aging(effects, dt=30*86400)  # Monthly updates
            
            # Update tissue composition
            self.update_tissue_composition(effects, dt=30*86400)
            
            # Update degeneration grade
            self.update_degeneration_grade()
            
    def simulate_injury_response(self, injury_params):
        """
        Simulate response to injury
        
        Parameters:
        -----------
        injury_params : dict
            Parameters describing the injury:
            - location: coordinates of injury center
            - severity: injury severity (0-1)
            - type: 'mechanical', 'chemical', or 'combined'
        """
        # Initialize injury response
        injury_location = injury_params['location']
        severity = injury_params['severity']
        
        # Calculate distance from injury site
        distances = np.sqrt(np.sum((self.nodes - injury_location)**2, axis=1))
        injury_field = severity * np.exp(-distances / self.params.INJURY_SPREAD)
        
        # Immediate effects
        if injury_params['type'] in ['mechanical', 'combined']:
            self.state['mechanical_damage'] += injury_field
            self.state['collagen_integrity'] *= (1 - 0.5 * injury_field)
            
        if injury_params['type'] in ['chemical', 'combined']:
            self.state['cytokines'] += injury_field * self.params.INJURY_CYTOKINE_FACTOR
            self.state['mmps'] += injury_field * self.params.INJURY_MMP_FACTOR
            
        # Simulate acute response phase
        acute_duration = 14  # days
        for day in range(acute_duration):
            self.update_injury_response(injury_field, dt=86400)
            
    def simulate_metabolic_stress(self, stress_params):
        """
        Simulate effects of metabolic stress
        
        Parameters:
        -----------
        stress_params : dict
            Parameters describing the metabolic stress:
            - oxygen_level: reduced oxygen level
            - glucose_level: reduced glucose level
            - ph_level: altered pH level
            - duration: stress duration in days
        """
        duration = stress_params['duration']
        
        for day in range(duration):
            # Update metabolic state
            self.update_metabolic_stress(stress_params, dt=86400)
            
            # Update cellular response
            self.update_cellular_metabolism(stress_params, dt=86400)
            
            # Update tissue response
            self.update_tissue_metabolism(stress_params, dt=86400)
            
    def update_cellular_aging(self, age_effects, dt):
        """Update cellular state during aging"""
        # Update senescent cell population
        senescence_increase = (age_effects['senescence_rate'] * 
                             (1 - self.state['senescent_cells']) * dt)
        self.state['senescent_cells'] += senescence_increase
        
        # Cell death
        cell_death = age_effects['cell_death_rate'] * dt
        self.state['senescent_cells'] *= (1 - cell_death)
        
        # Stem cell depletion
        self.state['stem_cells'] *= (1 - self.params.STEM_CELL_DEPLETION_RATE * dt)
        
    def update_tissue_composition(self, age_effects, dt):
        """Update tissue composition during aging"""
        # GAG loss
        gag_synthesis = age_effects['gag_synthesis'] * (1 - self.state['senescent_cells'])
        gag_degradation = self.params.BASE_GAG_DEGRADATION * (1 + self.state['mmps'])
        
        self.state['gag_content'] += (gag_synthesis - gag_degradation * 
                                     self.state['gag_content']) * dt
        
        # Collagen degradation
        collagen_damage = (self.params.BASE_COLLAGEN_DAMAGE * 
                          (1 + self.state['mechanical_damage']) * 
                          (1 + self.state['mmps']) * dt)
        
        self.state['collagen_integrity'] *= (1 - collagen_damage)
        
        # Water content changes
        target_water = 0.8 - 0.3 * (1 - self.state['gag_content'])
        water_adjustment = (target_water - self.state['water_content']) * 0.1
        self.state['water_content'] += water_adjustment * dt
        
    def update_injury_response(self, injury_field, dt):
        """Update tissue state during injury response"""
        # Inflammatory response
        inflammation = (self.state['cytokines'] * 
                      (1 + self.params.INJURY_AMPLIFICATION * injury_field))
        
        # MMP production
        mmp_production = (self.params.BASE_MMP_PRODUCTION * inflammation * 
                         (1 + self.state['senescent_cells']))
        mmp_degradation = self.params.MMP_DEGRADATION_RATE * self.state['mmps']
        
        self.state['mmps'] += (mmp_production - mmp_degradation) * dt
        
        # Growth factor response
        gf_production = self.params.INJURY_GF_PRODUCTION * (1 - inflammation)
        gf_degradation = self.params.GF_DEGRADATION_RATE * self.state['growth_factors']
        
        self.state['growth_factors'] += (gf_production - gf_degradation) * dt
        
        # Tissue healing
        healing_rate = (self.params.BASE_HEALING_RATE * 
                       self.state['growth_factors'] * 
                       (1 - inflammation))
        
        self.state['mechanical_damage'] *= (1 - healing_rate * dt)
        
    def update_metabolic_stress(self, stress_params, dt):
        """Update tissue state under metabolic stress"""
        # Calculate stress intensity
        o2_stress = max(0, 1 - stress_params['oxygen_level']/self.params.NORMAL_O2)
        glucose_stress = max(0, 1 - stress_params['glucose_level']/self.params.NORMAL_GLUCOSE)
        ph_stress = abs(stress_params['ph_level'] - 7.0)
        
        total_stress = (o2_stress + glucose_stress + ph_stress) / 3
        
        # Update oxidative stress
        ros_production = self.params.BASE_ROS_PRODUCTION * (1 + total_stress)
        self.state['oxidative_stress'] += ros_production * dt
        
        # Cell response to stress
        stress_induced_senescence = (self.params.STRESS_SENESCENCE_RATE * 
                                   total_stress * 
                                   (1 - self.state['senescent_cells']))
        
        self.state['senescent_cells'] += stress_induced_senescence * dt
        
    def update_degeneration_grade(self):
        """Update Pfirrmann grade based on tissue state"""
        # Calculate degeneration score
        score = (
            0.3 * (1 - self.state['water_content']) +
            0.3 * (1 - self.state['gag_content']) +
            0.2 * (1 - self.state['collagen_integrity']) +
            0.2 * self.state['senescent_cells']
        )
        
        # Map score to Pfirrmann grade (1-5)
        self.state['pfirrmann_grade'] = np.clip(1 + 4 * score, 1, 5)
        
    def get_degeneration_state(self):
        """Return current degeneration state"""
        return self.state.copy()

    def export_results(self, filename):
        """Export degeneration results to file"""
        results = {
            'state': self.state,
            'parameters': {
                name: getattr(self.params, name)
                for name in dir(self.params)
                if not name.startswith('_')
            }
        }
        
        np.savez(filename, **results)