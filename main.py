# main.py

import numpy as np
from datetime import datetime
import os
from ivd_model.modules.mechanical import MechanicalModel
from ivd_model.modules.biochemical import BiochemicalModel
from ivd_model.modules.cellular import CellularModel
from ivd_model.modules.spatial import SpatialModel
from ivd_model.modules.treatment import TreatmentModel
from ivd_model.validation.experimental_data import ExperimentalData
from ivd_model.validation.validation_metrics import ValidationMetrics
from ivd_model.analysis.sensitivity import SensitivityAnalysis
from ivd_model.analysis.uncertainty import UncertaintyAnalysis
from ivd_model.visualization.plotters import ResultsPlotter
from ivd_model.visualization.animations import AnimationGenerator
from scipy.interpolate import interp1d

class IVDSimulation:
    """
    Main simulation controller for IVD degeneration model
    """
    def __init__(self, simulation_name=None):
        # Create simulation name with timestamp
        if simulation_name is None:
            simulation_name = f"simulation_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        # Create results directory
        self.results_dir = f"results/{simulation_name}"
        os.makedirs(self.results_dir, exist_ok=True)
        
        # Initialize models
        n_points = 50   
        self.mechanical_model = MechanicalModel()
        self.biochemical_model = BiochemicalModel(spatial_discretization=n_points)
        self.cellular_model = CellularModel()
        self.spatial_model = SpatialModel(mesh_points=n_points)
        self.treatment_model = TreatmentModel()
        
        # Initialize analysis tools
        self.experimental_data = ExperimentalData()
        self.validation_metrics = ValidationMetrics()
        self.plotter = ResultsPlotter(save_dir=self.results_dir)
        self.animator = AnimationGenerator(save_dir=self.results_dir)
        
        # Initialize results storage
        self.results = {
            'time': [],
            'mechanical': [],
            'biochemical': [],
            'cellular': [],
            'spatial': []
        }

    def setup_sensitivity_analysis(self):
        """Setup sensitivity analysis parameters"""
        problem = {
            'num_vars': 6,
            'names': ['E_AF', 'E_NP', 'k_deg', 'D_O2', 'cell_density', 'applied_load'],
            'bounds': [
                [1.0, 4.0],    # E_AF (MPa)
                [0.5, 2.0],    # E_NP (MPa)
                [0.001, 0.01], # k_deg (1/day)
                [1e-3, 2e-3],  # D_O2 (mm²/s)
                [2e6, 6e6],    # cell_density (cells/mm³)
                [0.3, 1.5]     # applied_load (MPa)
            ]
        }
        return SensitivityAnalysis(problem)

    def setup_uncertainty_analysis(self):
        """Setup uncertainty analysis parameters"""
        uncertainty = UncertaintyAnalysis()
        
        # Add parameter distributions
        uncertainty.add_parameter_distribution('E_AF', 'normal', 
            {'mean': 2.5, 'std': 0.25})
        uncertainty.add_parameter_distribution('E_NP', 'normal',
            {'mean': 1.0, 'std': 0.1})
        uncertainty.add_parameter_distribution('k_deg', 'uniform',
            {'low': 0.001, 'high': 0.01})
            
        return uncertainty

    def run_simulation(self, days=365, applied_load=0.5, treatment=None):
        """Run main simulation"""
        dt = 1.0  # 1 day timestep
        time_points = np.arange(0, days+dt, dt)
        
        for t in time_points:
            # Update mechanical state
            mech_state = self.mechanical_model.solve_mechanics(applied_load)
            
            # Update biochemical state
            biochem_state = self.biochemical_model.update_concentrations(
                mech_state['height'],
                self.cellular_model.state['cell_density'],
                dt
            )
            
            # Update cellular state
            cell_state = self.cellular_model.update_state(
                biochem_state['oxygen'],
                biochem_state['glucose'],
                biochem_state['ph'],
                mech_state['stress'][0,2],  # Use axial stress
                dt
            )
            
            # Apply treatment if specified
            if treatment and t >= treatment['start_time']:
                treatment_effects = self.treatment_model.calculate_treatment_effect(
                    cell_state,
                    dt
                )
                cell_state.update(treatment_effects)
            
            # Update spatial distribution
            spatial_state = self.spatial_model.solve_mechanics(applied_load)
            
            # Store results
            self.results['time'].append(t)
            self.results['mechanical'].append(mech_state)
            self.results['biochemical'].append(biochem_state)
            self.results['cellular'].append(cell_state)
            self.results['spatial'].append(spatial_state)
            
            # Create progress visualization
            if int(t) % 30 == 0:  # Monthly updates
                self.plot_current_state(t)

    def validate_results(self):
        """Validate simulation results against experimental data"""
        # Prepare predictions
        predictions = {
            'height': [state['height'] for state in self.results['mechanical']],
            'water_content': [state['water_content'] for state in self.results['biochemical']],
            'cell_viability': [state['cell_viability'] for state in self.results['cellular']]
        }
        
        # Get experimental data
        experimental = {
            'height': self.experimental_data.get_validation_data('height'),
            'water_content': self.experimental_data.get_validation_data('water_content'),
            'cell_viability': self.experimental_data.get_validation_data('cell_viability')
        }
        
        # Calculate validation metrics
        validation_results = self.validation_metrics.validate_results(
            predictions,
            experimental
        )
        
        # Plot validation comparisons
        for metric in predictions.keys():
            self.plotter.plot_validation_comparison(
                predictions[metric],
                experimental[metric],
                f'{metric.replace("_", " ").title()} Validation',
                f'validation_{metric}.png'
            )
            
        return validation_results

    def run_sensitivity_analysis(self):
        """Perform sensitivity analysis"""
        sensitivity = self.setup_sensitivity_analysis()
        samples = sensitivity.generate_samples()
        
        # Run simulations for each sample
        outputs = []
        for sample in samples:
            # Set parameters
            self.mechanical_model.params.E_AF = sample[0]
            self.mechanical_model.params.E_NP = sample[1]
            self.cellular_model.params.k_DEG_COL = sample[2]
            self.biochemical_model.params.D_O2 = sample[3]
            self.cellular_model.params.CELL_DENSITY = sample[4]
            
            # Run simulation
            self.run_simulation(days=365, applied_load=sample[5])
            
            # Get output metric (e.g., final height loss)
            final_height = self.results['mechanical'][-1]['height']
            outputs.append(final_height)
            
        # Analyze sensitivity
        sensitivity_results = sensitivity.analyze_sensitivity(np.array(outputs))
        
        # Plot sensitivity results
        self.plotter.plot_sensitivity_analysis(
            sensitivity_results,
            'sensitivity_analysis.png'
        )
        
        return sensitivity_results

    def plot_current_state(self, time):
        """Plot current state of the simulation"""
        # Get current states
        mech_state = self.results['mechanical'][-1]
        biochem_state = self.results['biochemical'][-1]
        cell_state = self.results['cellular'][-1]
        
        # Create plots
        self.plotter.plot_time_series(
            self.results['time'],
            [state['height'] for state in self.results['mechanical']],
            'Disc Height Over Time',
            'Height (mm)',
            f'height_t{int(time)}.png'
        )
        
        # Get spatial data in correct format
        nodes = self.spatial_model.nodes
        oxygen = biochem_state['oxygen']
        
        # Ensure oxygen values match the number of nodes
        if len(oxygen) != len(nodes):
            # Interpolate oxygen to mesh points
            from scipy.interpolate import interp1d
            r = np.linspace(0, nodes[:,0].max(), len(oxygen))
            f = interp1d(r, oxygen, kind='cubic', fill_value='extrapolate')
            oxygen = f(nodes[:,0])
        
        self.plotter.plot_spatial_distribution(
            nodes[:,0],  # x coordinates
            nodes[:,1],  # y coordinates
            oxygen,      # interpolated oxygen values
            'Oxygen Distribution',
            f'oxygen_distribution_t{int(time)}.png'
        )

if __name__ == "__main__":
    # Create simulation instance
    sim = IVDSimulation()
    
    # Run baseline simulation
    print("Running baseline simulation...")
    sim.run_simulation(days=365, applied_load=0.5)
    
    # Validate results
    print("Validating results...")
    validation_results = sim.validate_results()
    
    # Run sensitivity analysis
    print("Performing sensitivity analysis...")
    sensitivity_results = sim.run_sensitivity_analysis()
    
    # Run treatment scenario
    print("Running treatment simulation...")
    treatment_params = {
        'start_time': 180,  # Start treatment at 6 months
        'type': 'cell_therapy',
        'parameters': {
            'cell_count': 1e6,
            'growth_factors': {
                'tgf_beta': 10.0,
                'igf1': 100.0
            }
        }
    }
    
    sim.run_simulation(days=365, applied_load=0.5, treatment=treatment_params)

    print("Creating animations...")
    sim.create_animations()
    
    print("Simulation completed successfully!")
    
    print("Simulation completed successfully!")