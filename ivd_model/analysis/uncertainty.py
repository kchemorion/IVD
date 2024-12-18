# analysis/uncertainty.py

import numpy as np
from scipy.stats import norm, uniform, multivariate_normal

class UncertaintyAnalysis:
    """
    Uncertainty quantification and propagation
    """
    
    def __init__(self):
        self.parameter_distributions = {}
        
    def add_parameter_distribution(self, name, dist_type, params):
        """Add parameter distribution"""
        self.parameter_distributions[name] = {
            'type': dist_type,
            'params': params
        }
        
    def sample_parameters(self, N=1000):
        """Generate parameter samples"""
        samples = {}
        for name, dist in self.parameter_distributions.items():
            if dist['type'] == 'normal':
                samples[name] = norm.rvs(loc=dist['params']['mean'],
                                       scale=dist['params']['std'],
                                       size=N)
            elif dist['type'] == 'uniform':
                samples[name] = uniform.rvs(loc=dist['params']['low'],
                                          scale=dist['params']['high']-dist['params']['low'],
                                          size=N)
        return samples
    
    def calculate_uncertainty_bounds(self, results, confidence=0.95):
        """Calculate confidence intervals"""
        bounds = {}
        for metric, values in results.items():
            values = np.array(values)
            mean = np.mean(values, axis=0)
            std = np.std(values, axis=0)
            
            z = norm.ppf((1 + confidence) / 2)
            lower = mean - z * std
            upper = mean + z * std
            
            bounds[metric] = {
                'mean': mean,
                'lower': lower,
                'upper': upper,
                'std': std
            }
        return bounds