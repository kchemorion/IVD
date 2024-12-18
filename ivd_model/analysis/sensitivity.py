# analysis/sensitivity.py

import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol

class SensitivityAnalysis:
    """
    Sensitivity analysis using Sobol method
    """
    
    def __init__(self, problem_dict):
        """
        Initialize with problem dictionary defining parameters
        Example:
        {
            'num_vars': 3,
            'names': ['E_AF', 'E_NP', 'k_deg'],
            'bounds': [[1.0, 4.0], [0.5, 2.0], [0.001, 0.01]]
        }
        """
        self.problem = problem_dict
        
    def generate_samples(self, N=1000):
        """Generate parameter samples"""
        return saltelli.sample(self.problem, N)
    
    def analyze_sensitivity(self, Y):
        """Analyze sensitivity using Sobol indices"""
        Si = sobol.analyze(self.problem, Y)
        return {
            'S1': Si['S1'],
            'ST': Si['ST'],
            'S2': Si['S2'],
            'names': self.problem['names']
        }
    
    def rank_parameters(self, sensitivity_indices):
        """Rank parameters by importance"""
        total_effects = sensitivity_indices['ST']
        names = sensitivity_indices['names']
        
        # Sort by total effect
        ranking = sorted(zip(names, total_effects), 
                       key=lambda x: x[1], reverse=True)
        return ranking