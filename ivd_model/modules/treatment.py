# modules/treatment.py

import numpy as np
from ..config.parameters import TreatmentParams as TP

class TreatmentModel:
    """
    Treatment intervention model for IVD regeneration
    Based on:
    1. Sakai & Andersson (2015) - DOI: 10.1038/nrrheum.2014.198
    2. Richardson et al. (2016) - DOI: 10.1038/nmat4614
    """
    
    def __init__(self):
        self.params = TP
        self.reset_state()
        
    def reset_state(self):
        self.state = {
            'msc_count': 0,
            'growth_factors': {
                'tgf_beta': 0.0,
                'igf1': 0.0,
                'bmp2': 0.0
            },
            'scaffold_properties': {
                'porosity': 0.0,
                'stiffness': 0.0,
                'degradation_rate': 0.0
            },
            'treatment_active': False
        }
        
    def apply_cell_therapy(self, cell_count):
        """Apply MSC-based cell therapy"""
        survival_rate = self.params.MSC_SURVIVAL
        diff_rate = self.params.MSC_DIFFERENTIATION
        
        self.state['msc_count'] = cell_count
        self.state['treatment_active'] = True
        
        return {
            'surviving_cells': cell_count * survival_rate,
            'differentiated_cells': cell_count * survival_rate * diff_rate
        }
    
    def apply_growth_factors(self, factors):
        """Apply growth factor therapy"""
        self.state['growth_factors'].update(factors)
        self.state['treatment_active'] = True
        
    def apply_scaffold(self, properties):
        """Apply biomaterial scaffold"""
        self.state['scaffold_properties'].update(properties)
        self.state['treatment_active'] = True
        
    def calculate_treatment_effect(self, current_state, dt):
        """Calculate effect of current treatment"""
        if not self.state['treatment_active']:
            return current_state
            
        effects = {
            'cell_density': current_state['cell_density'],
            'collagen': current_state['collagen'],
            'proteoglycan': current_state['proteoglycan'],
            'water_content': current_state['water_content']
        }
        
        # Cell therapy effects
        if self.state['msc_count'] > 0:
            surviving_cells = self.state['msc_count'] * np.exp(-0.05 * dt)  # Cell death rate
            effects['cell_density'] += surviving_cells
            self.state['msc_count'] = surviving_cells
            
        # Growth factor effects
        if self.state['growth_factors']['tgf_beta'] > 0:
            tgf_effect = self.params.TGF_BETA_FACTOR
            effects['collagen'] *= tgf_effect
            effects['proteoglycan'] *= tgf_effect
            
        if self.state['growth_factors']['igf1'] > 0:
            igf_effect = self.params.IGF1_FACTOR
            effects['proteoglycan'] *= igf_effect
            
        # Scaffold effects
        if any(self.state['scaffold_properties'].values()):
            scaffold = self.state['scaffold_properties']
            mechanical_support = scaffold['stiffness'] / 100.0  # Normalize to 0-1
            effects['water_content'] *= (1 + scaffold['porosity']/2)
            
        return effects