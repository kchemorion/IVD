# validation/experimental_data.py

import numpy as np
import pandas as pd
from ..config.validation_data import *

class ExperimentalData:
    """Handle experimental validation data"""
    
    def __init__(self):
        self.height_data = pd.DataFrame(HEIGHT_DATA)
        self.water_content_data = pd.DataFrame(WATER_CONTENT_DATA)
        self.cell_viability_data = pd.DataFrame(CELL_VIABILITY_DATA)
        self.proteoglycan_data = pd.DataFrame(PROTEOGLYCAN_DATA)
        self.mechanical_data = MECHANICAL_VALIDATION
        self.treatment_data = TREATMENT_VALIDATION
        
    def get_validation_data(self, data_type):
        """Get validation data for specific metric"""
        if data_type == 'height':
            return {
                'time_points': HEIGHT_DATA['time_points'],
                'values': HEIGHT_DATA['heights'],
                'std_dev': HEIGHT_DATA['std_dev']
            }
        elif data_type == 'water_content':
            return {
                'time_points': WATER_CONTENT_DATA['time_points'],
                'values': WATER_CONTENT_DATA['content'],
                'std_dev': WATER_CONTENT_DATA['std_dev']
            }
        elif data_type == 'cell_viability':
            return {
                'time_points': CELL_VIABILITY_DATA['time_points'],
                'values': CELL_VIABILITY_DATA['viability'],
                'std_dev': CELL_VIABILITY_DATA['std_dev']
            }
        else:
            raise ValueError(f"Unknown data type: {data_type}")