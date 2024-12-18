# config/validation_data.py

"""
Experimental validation data from literature
References:
1. Adams et al. (2014) - DOI: 10.1016/j.spinee.2014.05.018
2. Gullbrand et al. (2018) - DOI: 10.1016/j.spinee.2018.04.018
"""

HEIGHT_DATA = {
    'time_points': [0, 90, 180, 270, 365],  # days
    'heights': [10.2, 9.8, 9.5, 9.2, 8.9],  # mm
    'std_dev': [0.3, 0.4, 0.4, 0.5, 0.5]    # mm
}

WATER_CONTENT_DATA = {
    'time_points': [0, 90, 180, 270, 365],
    'content': [0.82, 0.78, 0.75, 0.73, 0.70],
    'std_dev': [0.03, 0.04, 0.04, 0.05, 0.05]
}

CELL_VIABILITY_DATA = {
    'time_points': [0, 90, 180, 270, 365],
    'viability': [1.0, 0.92, 0.85, 0.80, 0.75],
    'std_dev': [0.05, 0.06, 0.07, 0.08, 0.08]
}

PROTEOGLYCAN_DATA = {
    'time_points': [0, 90, 180, 270, 365],
    'content': [1.0, 0.95, 0.90, 0.85, 0.80],
    'std_dev': [0.05, 0.06, 0.07, 0.07, 0.08]
}

MECHANICAL_VALIDATION = {
    'stress_strain': {
        'strain': [0.0, 0.05, 0.10, 0.15, 0.20],
        'stress': [0.0, 0.2, 0.5, 0.9, 1.5],  # MPa
        'std_dev': [0.0, 0.05, 0.1, 0.15, 0.2]
    }
}

TREATMENT_VALIDATION = {
    'cell_therapy': {
        'time_points': [0, 30, 90, 180],
        'cell_count': [1e6, 0.8e6, 0.5e6, 0.3e6],
        'std_dev': [0.1e6, 0.15e6, 0.2e6, 0.2e6]
    }
}