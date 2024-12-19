"""
Parameters for IVD degeneration model
Based on experimental data from literature
"""

class DegenerationParams:
    # Age-related parameters
    BASE_CELL_DEATH_RATE = 0.001      # 1/day
    AGE_DEATH_FACTOR = 0.02           # Increase per year above 20
    BASE_GAG_SYNTHESIS = 1e-5         # g/ml/day
    AGE_GAG_DECAY = 0.03              # Decay rate with age
    BASE_SENESCENCE_RATE = 0.0005     # 1/day
    AGE_SENESCENCE_FACTOR = 0.03      # Increase per year above 20
    STEM_CELL_DEPLETION_RATE = 0.0002 # 1/day
    
    # Tissue composition parameters
    BASE_GAG_DEGRADATION = 0.001      # 1/day
    BASE_COLLAGEN_DAMAGE = 0.0005     # 1/day
    NORMAL_WATER_CONTENT = 0.8        # fraction
    MIN_WATER_CONTENT = 0.5           # fraction
    
    # Injury response parameters
    INJURY_SPREAD = 0.005             # meters
    INJURY_CYTOKINE_FACTOR = 10.0     # relative to baseline
    INJURY_MMP_FACTOR = 5.0           # relative to baseline
    INJURY_AMPLIFICATION = 2.0        # inflammation amplification
    BASE_MMP_PRODUCTION = 0.1         # concentration/day
    MMP_DEGRADATION_RATE = 0.5        # 1/day
    INJURY_GF_PRODUCTION = 0.2        # concentration/day
    GF_DEGRADATION_RATE = 0.3         # 1/day
    BASE_HEALING_RATE = 0.1           # 1/day
    
    # Metabolic stress parameters
    NORMAL_O2 = 5.0                   # kPa
    NORMAL_GLUCOSE = 4.0              # mM
    BASE_ROS_PRODUCTION = 0.01        # concentration/day
    STRESS_SENESCENCE_RATE = 0.01     # 1/day
    
    # Recovery parameters
    MAX_REGENERATION_RATE = 0.05      # 1/day
    REGENERATION_THRESHOLD = 0.3      # minimum growth factors needed
    
    # Critical thresholds
    CRITICAL_GAG_CONTENT = 0.3        # fraction of normal
    CRITICAL_COLLAGEN_INTEGRITY = 0.4  # fraction of normal
    CRITICAL_CELL_DENSITY = 0.2       # fraction of normal