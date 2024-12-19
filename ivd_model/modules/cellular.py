import numpy as np
from scipy.integrate import solve_ivp
from ..config.parameters import CellularParams as CP

class CellularModel:
    def update_state(self, dt, mechanical_input, biochemical_input):
        """
        Update cellular state for one time step
        
        Parameters:
        -----------
        dt : float
            Time step size
        mechanical_input : dict
            Current mechanical state (stress, strain)
        biochemical_input : dict
            Current biochemical state (oxygen, glucose, pH)
        """
        # Update stimuli
        self.calculate_mechanical_stimulus(
            mechanical_input['stress'],
            mechanical_input['strain']
        )
        
        self.calculate_metabolic_stress(
            biochemical_input['oxygen'],
            biochemical_input['glucose'],
            biochemical_input['ph']
        )
        
        # Update cell populations
        self.update_cell_populations(dt)
        
        # Update ECM synthesis and degradation
        self.update_ecm(dt)
        
        # Update inflammatory mediators
        self.update_inflammation(dt)
        
        # Update geometry visualization
        self.update_geometry()
        
    def update_cell_populations(self, dt):
        """Update cell populations based on current conditions"""
        # Get current states
        mech_stim = self.state['mechanical_stimulus']
        metab_stress = self.state['metabolic_stress']
        
        # Base rates
        death_rate = self.transition_rates['senescent_to_death']
        senescence_rate = self.transition_rates['healthy_to_senescent']
        prolif_rate = self.transition_rates['proliferation']
        
        # Modify rates based on environmental conditions
        death_rate *= (1 + metab_stress)
        senescence_rate *= (1 + 0.5*metab_stress + 0.5*np.abs(mech_stim))
        prolif_rate *= (1 - metab_stress) * np.clip(1 + mech_stim, 0, 2)
        
        # Update populations
        healthy = self.state['healthy_fraction']
        senescent = self.state['senescent_fraction']
        density = self.state['cell_density']
        
        # Population dynamics
        dH = (-senescence_rate * healthy + prolif_rate * healthy * 
              (1 - density/self.params.CELL_DENSITY_MAX))
        dS = (senescence_rate * healthy - death_rate * senescent)
        
        # Update fractions
        self.state['healthy_fraction'] += dH * dt
        self.state['senescent_fraction'] += dS * dt
        
        # Update total density
        total_change = (prolif_rate * healthy * 
                       (1 - density/self.params.CELL_DENSITY_MAX) - 
                       death_rate * senescent) * density
        self.state['cell_density'] += total_change * dt
        
        # Ensure fractions sum to 1 and density is positive
        total_fraction = self.state['healthy_fraction'] + self.state['senescent_fraction']
        self.state['healthy_fraction'] /= total_fraction
        self.state['senescent_fraction'] /= total_fraction
        np.clip(self.state['cell_density'], 0, self.params.CELL_DENSITY_MAX, 
                out=self.state['cell_density'])
        
    def update_ecm(self, dt):
        """Update ECM composition"""
        # Get cell populations
        healthy = self.state['healthy_fraction']
        senescent = self.state['senescent_fraction']
        density = self.state['cell_density']
        
        # Calculate synthesis rates
        collagen_synth = (
            healthy * self.phenotype_properties['healthy']['ecm_synthesis'] +
            senescent * self.phenotype_properties['senescent']['ecm_synthesis']
        ) * density
        
        # MMP production
        mmp_prod = (
            healthy * self.phenotype_properties['healthy']['mmp_production'] +
            senescent * self.phenotype_properties['senescent']['mmp_production']
        ) * density
        
        # Update MMPs and TIMPs
        self.state['mmp'] += (mmp_prod - self.params.MMP_DECAY * self.state['mmp']) * dt
        self.state['timp'] += (self.params.TIMP_PRODUCTION * density - 
                              self.params.TIMP_DECAY * self.state['timp']) * dt
        
        # Calculate degradation rates
        mmp_active = np.maximum(self.state['mmp'] - self.state['timp'], 0)
        degradation_rate = self.params.ECM_DEGRADATION * mmp_active
        
        # Update ECM components
        self.state['collagen'] += (collagen_synth - 
                                  degradation_rate * self.state['collagen']) * dt
        
        self.state['proteoglycan'] += (self.params.PG_SYNTHESIS * density - 
                                      degradation_rate * self.state['proteoglycan']) * dt
        
        self.state['crosslinks'] += (self.params.CROSSLINK_FORMATION * self.state['collagen'] - 
                                    self.params.CROSSLINK_DEGRADATION * self.state['crosslinks']) * dt
        
    def update_inflammation(self, dt):
        """Update inflammatory mediators"""
        # Get cell populations
        healthy = self.state['healthy_fraction']
        senescent = self.state['senescent_fraction']
        density = self.state['cell_density']
        
        # Calculate inflammatory mediator production
        il1_prod = (
            healthy * self.phenotype_properties['healthy']['inflammation'] +
            senescent * self.phenotype_properties['senescent']['inflammation']
        ) * density
        
        tnf_prod = il1_prod * self.params.TNF_IL1_RATIO
        
        # Update cytokine concentrations
        self.state['il1'] += (il1_prod - self.params.IL1_DECAY * self.state['il1']) * dt
        self.state['tnf_alpha'] += (tnf_prod - self.params.TNF_DECAY * self.state['tnf_alpha']) * dt
        
        # Update inflammatory state
        self.state['inflammatory_state'] = (
            self.state['il1'] / self.params.IL1_REFERENCE +
            self.state['tnf_alpha'] / self.params.TNF_REFERENCE
        ) / 2
        
    def calculate_ecm_properties(self):
        """Calculate ECM mechanical properties based on composition"""
        # Normalize ECM components
        collagen_norm = self.state['collagen'] / self.params.COLLAGEN_INITIAL
        pg_norm = self.state['proteoglycan'] / self.params.PG_INITIAL
        crosslinks_norm = self.state['crosslinks'] / self.params.CROSSLINKS_INITIAL
        
        # Calculate properties
        properties = {
            'matrix_stiffness': (
                self.params.BASE_STIFFNESS * 
                (collagen_norm * (1 + self.params.CROSSLINK_STIFFNESS_FACTOR * crosslinks_norm))
            ),
            'swelling_pressure': (
                self.params.BASE_SWELLING * pg_norm
            ),
            'permeability': (
                self.params.BASE_PERMEABILITY * 
                np.exp(-self.params.PERM_COLLAGEN_FACTOR * collagen_norm)
            )
        }
        
        return properties
        
    def calculate_metabolic_rates(self, biochemical_state):
        """Calculate metabolic consumption/production rates"""
        # Get cell density and viability
        density = self.state['cell_density']
        healthy = self.state['healthy_fraction']
        
        # Base metabolic rates
        o2_consumption = self.params.O2_CONSUMPTION_RATE * density * healthy
        glucose_consumption = self.params.GLUCOSE_CONSUMPTION_RATE * density * healthy
        lactate_production = self.params.LACTATE_PRODUCTION_RATE * density * healthy
        
        # Modify by metabolic state
        hypoxia_factor = np.clip(
            biochemical_state['oxygen'] / self.params.O2_CRITICAL, 
            0, 1
        )
        
        # Return consumption rates
        return {
            'oxygen': -o2_consumption * hypoxia_factor,
            'glucose': -glucose_consumption * hypoxia_factor,
            'lactate': lactate_production * (2 - hypoxia_factor)
        }
        
    def update_geometry(self):
        """Update geometry model with current solution"""
        for field in self.state:
            self.geometry.interpolate_to_vertices(
                field,
                self.state[field],
                self.nodes
            )
