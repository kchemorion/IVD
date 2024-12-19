"""
Experimental validation data from literature sources
References:
1. Urban et al. (2004) - Nutrition of the intervertebral disc
2. Galbusera et al. (2011) - Ageing and degenerative changes of the IVD
3. Johannessen et al. (2006) - Disc mechanics during daily living
4. Iatridis et al. (2013) - Mechanical behavior and structure of the disc
"""

import numpy as np
import pandas as pd

class ExperimentalData:
    def __init__(self):
        # Load mechanical validation data
        self.mechanical_data = {
            'compression': {
                'strain': np.array([0, 0.05, 0.10, 0.15, 0.20]),
                'stress': np.array([0, 0.15, 0.45, 0.90, 1.80]),  # MPa
                'source': 'Johannessen et al. 2006'
            },
            'creep': {
                'time': np.array([0, 1, 2, 4, 8, 16, 24]),  # hours
                'displacement': np.array([0, 0.8, 1.2, 1.5, 1.7, 1.85, 1.9]),  # mm
                'source': 'Iatridis et al. 2013'
            }
        }

        # Biochemical concentration profiles
        self.biochemical_data = {
            'oxygen': {
                'radius': np.linspace(0, 25, 10),  # mm from center
                'concentration': np.array([1.5, 2.0, 2.8, 3.5, 4.2, 5.0, 5.8, 6.2, 6.8, 7.0]),  # kPa
                'source': 'Urban et al. 2004'
            },
            'glucose': {
                'radius': np.linspace(0, 25, 10),  # mm from center
                'concentration': np.array([2.0, 2.5, 3.0, 3.3, 3.8, 4.0, 4.3, 4.5, 4.8, 5.0]),  # mM
                'source': 'Urban et al. 2004'
            }
        }

        # Cell viability data
        self.cellular_data = {
            'age_related': {
                'age': np.array([20, 30, 40, 50, 60, 70]),  # years
                'density': np.array([4e6, 3.8e6, 3.5e6, 3e6, 2.5e6, 2e6]),  # cells/ml
                'source': 'Galbusera et al. 2011'
            },
            'degeneration': {
                'grade': np.array([1, 2, 3, 4, 5]),  # Pfirrmann grade
                'viability': np.array([0.95, 0.85, 0.70, 0.50, 0.30]),  # fraction
                'source': 'Galbusera et al. 2011'
            }
        }

        # Treatment response data
        self.treatment_data = {
            'cell_therapy': {
                'time': np.array([0, 1, 2, 4, 8, 12, 24, 48]),  # weeks
                'cell_survival': np.array([1.0, 0.8, 0.65, 0.5, 0.35, 0.25, 0.15, 0.10]),
                'source': 'Literature compilation'
            },
            'growth_factors': {
                'time': np.array([0, 1, 2, 4, 7, 14, 21, 28]),  # days
                'concentration': np.array([1.0, 0.7, 0.5, 0.3, 0.15, 0.08, 0.04, 0.02]),
                'source': 'Literature compilation'
            }
        }

    def get_validation_data(self, data_type, subtype=None):
        """
        Retrieve validation data for specific aspect of the model
        
        Parameters:
        -----------
        data_type : str
            Type of data ('mechanical', 'biochemical', 'cellular', 'treatment')
        subtype : str, optional
            Specific subtype of data
            
        Returns:
        --------
        dict
            Dictionary containing validation data
        """
        if data_type == 'mechanical':
            return self.mechanical_data if subtype is None else self.mechanical_data[subtype]
        elif data_type == 'biochemical':
            return self.biochemical_data if subtype is None else self.biochemical_data[subtype]
        elif data_type == 'cellular':
            return self.cellular_data if subtype is None else self.cellular_data[subtype]
        elif data_type == 'treatment':
            return self.treatment_data if subtype is None else self.treatment_data[subtype]
        else:
            raise ValueError(f"Unknown data type: {data_type}")

    def calculate_rmse(self, predictions, data_type, subtype=None):
        """Calculate Root Mean Square Error between predictions and experimental data"""
        validation_data = self.get_validation_data(data_type, subtype)
        
        if isinstance(predictions, np.ndarray) and isinstance(validation_data, dict):
            if 'concentration' in validation_data:
                experimental = validation_data['concentration']
            elif 'stress' in validation_data:
                experimental = validation_data['stress']
            elif 'displacement' in validation_data:
                experimental = validation_data['displacement']
            else:
                raise ValueError("Incompatible data format")
                
            return np.sqrt(np.mean((predictions - experimental) ** 2))
        else:
            raise ValueError("Invalid prediction format")

    def validate_model(self, model_predictions):
        """
        Comprehensive model validation against all available experimental data
        
        Parameters:
        -----------
        model_predictions : dict
            Dictionary containing model predictions for different aspects
            
        Returns:
        --------
        dict
            Dictionary containing validation metrics for each aspect
        """
        validation_results = {}
        
        # Validate mechanical response
        if 'mechanical' in model_predictions:
            mech_rmse = self.calculate_rmse(
                model_predictions['mechanical']['stress'],
                'mechanical',
                'compression'
            )
            validation_results['mechanical_rmse'] = mech_rmse
        
        # Validate biochemical concentrations
        if 'biochemical' in model_predictions:
            for species in ['oxygen', 'glucose']:
                if species in model_predictions['biochemical']:
                    bio_rmse = self.calculate_rmse(
                        model_predictions['biochemical'][species],
                        'biochemical',
                        species
                    )
                    validation_results[f'{species}_rmse'] = bio_rmse
        
        # Validate cell viability
        if 'cellular' in model_predictions and 'viability' in model_predictions['cellular']:
            cell_rmse = self.calculate_rmse(
                model_predictions['cellular']['viability'],
                'cellular',
                'degeneration'
            )
            validation_results['cell_viability_rmse'] = cell_rmse
        
        return validation_results