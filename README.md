# Advanced Intervertebral Disc (IVD) Computational Model

## Overview
This repository contains a comprehensive computational model of the intervertebral disc that integrates mechanical, biochemical, and cellular processes. The model is designed to simulate various aspects of IVD behavior, degeneration, and therapeutic responses using real geometric data and validated mathematical models from literature.

## Key Features

### 1. Geometric Representation
- Uses real IVD geometry from medical imaging (STL format)
- Distinguishes between Annulus Fibrosus (AF) and Nucleus Pulposus (NP) regions
- Supports patient-specific geometry adaptation
- Implements proper mesh handling and refinement capabilities

### 2. Mechanical Model
References: [Adams & Roughley, 2006](doi:10.1097/01.brs.0000231761.73859.2c), [Iatridis et al., 2006](doi:10.1007/s00586-006-0109-9)
- Anisotropic biphasic material properties
- Fiber-reinforced constitutive model
- Region-specific mechanical properties
- Handles complex loading scenarios
- Includes poroelastic behavior
- Accounts for fluid flow and pressure distribution

### 3. Biochemical Transport
References: [Urban & Roberts, 2003](doi:10.1186/ar629), [Bibby et al., 2005](doi:10.1097/01.brs.0000171207.78927.d9)
- 3D reaction-diffusion equations
- Models key molecules:
  - Oxygen
  - Glucose
  - Lactate
  - Growth factors
  - Inflammatory mediators
- pH dynamics
- Coupled with mechanical deformation

### 4. Cellular Components
References: [Huang et al., 2014](doi:10.1038/nrrheum.2014.91), [Rodriguez et al., 2011](doi:10.1111/j.1474-9726.2011.00700.x)
- Cell viability and density dynamics
- Mechanotransduction
- Metabolic responses
- ECM synthesis and degradation
- Cell senescence
- Inflammatory responses

### 5. Degeneration Mechanisms
References: [Vergroesen et al., 2015](doi:10.1016/j.joca.2015.05.017), [Antoniou et al., 1996](doi:10.1097/00007632-199610010-00006)

#### Age-related Degeneration
- Progressive tissue changes
- Cell senescence
- ECM degradation
- Water content changes
- Mechanical property alterations

#### Injury-induced Degeneration
- Annular tears
- Endplate damage
- Mechanical trauma
- Inflammatory cascade
- Tissue remodeling

#### Metabolic Stress
- Nutrient deprivation
- pH changes
- Oxidative stress
- Cell death
- Matrix degradation

### 6. Treatment Modeling
- Cell therapy optimization
- Growth factor delivery
- Biomaterial applications
- Combined treatment strategies
- Therapeutic response prediction

### 7. Visualization Capabilities
- 3D visualization of results on real geometry
- Time-series analysis
- Cross-sectional views
- Interactive dashboard
- Export capabilities for various formats

## Implementation Details

### Directory Structure
```
ivd_model/
├── geometry/         # Geometry handling and mesh operations
├── modules/          # Core simulation modules
│   ├── mechanical.py
│   ├── biochemical.py
│   ├── cellular.py
│   ├── degeneration.py
│   └── treatment.py
├── config/          # Model parameters and configuration
├── visualization/   # Visualization tools
├── validation/     # Experimental validation data
├── studies/        # Case studies and analysis
└── tests/          # Testing suite
```

### Key Dependencies
- Python ≥ 3.8
- NumPy ≥ 1.20.0
- SciPy ≥ 1.7.0
- Pandas ≥ 1.3.0
- Meshio ≥ 5.0.0
- Plotly ≥ 5.0.0

## Validation

The model has been validated against experimental data from literature for:
1. Mechanical response under various loading conditions
2. Nutrient distribution patterns
3. Cell viability profiles
4. Degeneration progression
5. Treatment outcomes

### Key Validation Sources
- Mechanical properties: [Johannessen et al., 2006](doi:10.1097/01.brs.0000231761.73859.2c)
- Transport properties: [Urban et al., 2004](doi:10.1097/01.brs.0000231761.73859.2c)
- Cell behavior: [Bibby & Urban, 2004](doi:10.1097/01.brs.0000142434.61410.33)
- Degeneration: [Antoniou et al., 1996](doi:10.1097/00007632-199610010-00006)

## Usage

### Basic Usage
```python
from ivd_model import IVDModel

# Initialize model with geometry
model = IVDModel("path/to/geometry.stl")

# Run simulation
results = model.simulate(days=365)

# Visualize results
model.visualize(results)
```

### Running Case Studies
```python
from ivd_model.studies import DegenerationStudies

# Initialize studies
studies = DegenerationStudies(model)

# Run age-related degeneration study
results = studies.run_age_related_study(
    starting_age=20,
    duration=40  # years
)
```

## Contributing
1. Fork the repository
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Create a new Pull Request

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use this model in your research, please cite our work:
[Citation details to be added after publication]

## Authors
[Author names and affiliations]

## Acknowledgments
- Model development supported by [funding sources]
- Geometry data provided by [source]
- Validation data from multiple published sources (see references)

## References
[Complete list of references used in model development and validation]
