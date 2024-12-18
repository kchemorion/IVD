# config/parameters.py

"""
IVD Model Parameters
References:
1. Mechanical parameters: Malandrino et al. (2015) DOI: 10.1016/j.jbiomech.2015.02.002
2. Biochemical parameters: Neidlinger-Wilke et al. (2014) DOI: 10.1016/j.jbbm.2013.04.010
3. Cell metabolism: Huang & Gu (2008) DOI: 10.1016/j.jbiomech.2008.03.009
"""

class MechanicalParams:
    # Tissue properties
    E_AF = 2.5        # AF Young's modulus (MPa)
    E_NP = 1.0        # NP Young's modulus (MPa)
    v_AF = 0.325      # AF Poisson's ratio
    v_NP = 0.49       # NP Poisson's ratio
    k0_AF = 0.0025    # AF initial permeability (mm^4/Ns)
    k0_NP = 0.001     # NP initial permeability (mm^4/Ns)
    
    # Geometry
    DISC_HEIGHT = 10.0  # Initial disc height (mm)
    AF_RATIO = 0.7      # AF outer/inner radius ratio
    
    # Osmotic properties
    c_FCD = -0.2      # Fixed charge density (mEq/mm^3)
    phi0 = 0.8        # Initial water content

class BiochemicalParams:
    # Diffusion coefficients (mm^2/s)
    D_O2 = 1.5e-3     # Oxygen
    D_GLU = 1.0e-3    # Glucose
    D_LAC = 1.0e-3    # Lactate
    
    # Consumption/production rates
    Q_O2 = 1.8e-8     # Oxygen consumption (mol/cell/day)
    Q_GLU = 5.6e-8    # Glucose consumption (mol/cell/day)
    Q_LAC = 1.12e-7   # Lactate production (mol/cell/day)
    
    # Concentration bounds
    O2_CRIT = 1.0     # Critical oxygen tension (kPa)
    GLU_CRIT = 0.5    # Critical glucose concentration (mmol/L)
    LAC_MAX = 10.0    # Maximum lactate concentration (mmol/L)

class CellularParams:
    # Cell properties
    CELL_DENSITY = 4e6  # Cells per mm^3
    CELL_RADIUS = 10.0  # Cell radius (μm)
    
    # Synthesis rates (/day)
    k_COL = 2.3e-3     # Collagen
    k_PG = 2.7e-3      # Proteoglycan
    
    # Degradation rates (/day)
    k_DEG_COL = 1.2e-3  # Collagen
    k_DEG_PG = 1.5e-3   # Proteoglycan

class TreatmentParams:
    # Growth factor effects
    TGF_BETA_FACTOR = 2.0    # TGF-β effect on synthesis
    IGF1_FACTOR = 1.5        # IGF-1 effect on synthesis
    
    # Cell therapy parameters
    MSC_SURVIVAL = 0.3       # MSC survival rate
    MSC_DIFFERENTIATION = 0.5 # MSC differentiation rate

class SimulationParams:
    # Time settings
    TIME_STEPS_PER_DAY = 24
    TOTAL_DAYS = 365
    DT = 1.0/TIME_STEPS_PER_DAY
    
    # Numerical parameters
    CONV_TOL = 1e-6
    MAX_ITERATIONS = 100
    
    # Output settings
    OUTPUT_FREQUENCY = 1  # days

# Validation thresholds
VALIDATION_THRESHOLDS = {
    'height_rmse': 0.5,      # mm
    'stress_rmse': 0.2,      # MPa
    'water_content_rmse': 0.05,
    'cell_viability_rmse': 0.1
}