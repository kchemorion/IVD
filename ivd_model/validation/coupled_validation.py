import numpy as np
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

class CoupledValidation:
    """
    Validation tools for coupled IVD model behavior.
    Compares simulation results with experimental data and literature values.
    """
    
    def __init__(self, coupling_model):
        self.coupling = coupling_model
        self.experimental_data = {}
        self.validation_metrics = {}
        self.simulation_results = None
        self.load_experimental_data()
        
    def load_experimental_data(self):
        """Load experimental validation data"""
        # This would typically load from files or a database
        # For now, using hardcoded literature values from papers
        self.experimental_data = {
            'mechanical': {
                'compression': {
                    'strain': np.array([0, 0.05, 0.10, 0.15]),
                    'stress': np.array([0, 0.2, 0.5, 1.0])  # MPa
                },
                'creep': {
                    'time': np.array([0, 1, 2, 4, 8, 24]),  # hours
                    'displacement': np.array([0, 0.5, 0.8, 1.1, 1.3, 1.5])  # mm
                }
            },
            'biochemical': {
                'oxygen': {
                    'radial_position': np.linspace(0, 1, 10),
                    'concentration': np.array([5.8, 5.2, 4.5, 3.8, 3.1, 2.5, 2.0, 1.5, 1.1, 0.8])
                },
                'glucose': {
                    'radial_position': np.linspace(0, 1, 10),
                    'concentration': np.array([4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.2, 0.9, 0.6, 0.4])
                }
            },
            'cellular': {
                'viability': {
                    'oxygen_level': np.array([0.1, 0.5, 1.0, 2.0, 3.0, 5.0]),
                    'viability': np.array([0.1, 0.3, 0.6, 0.8, 0.9, 1.0])
                }
            }
        }
        
    def validate_mechanical_response(self, simulation_results):
        """Validate mechanical behavior"""
        self.simulation_results = simulation_results
        exp_data = self.experimental_data['mechanical']
        
        # Validate compression response
        sim_stress = np.interp(
            exp_data['compression']['strain'],
            simulation_results['strain'],
            simulation_results['stress']
        )
        
        stress_error = np.mean(
            np.abs(sim_stress - exp_data['compression']['stress']) /
            exp_data['compression']['stress']
        )
        
        stress_correlation, _ = pearsonr(sim_stress, 
                                       exp_data['compression']['stress'])
        
        # Validate creep response
        sim_displacement = np.interp(
            exp_data['creep']['time'],
            simulation_results['time'],
            simulation_results['displacement']
        )
        
        displacement_error = np.mean(
            np.abs(sim_displacement - exp_data['creep']['displacement']) /
            exp_data['creep']['displacement']
        )
        
        self.validation_metrics['mechanical'] = {
            'stress_error': stress_error,
            'stress_correlation': stress_correlation,
            'displacement_error': displacement_error
        }
        
        return self.validation_metrics['mechanical']
        
    def validate_biochemical_distribution(self, simulation_results):
        """Validate biochemical distributions"""
        self.simulation_results = simulation_results
        exp_data = self.experimental_data['biochemical']
        metrics = {}
        
        for species in ['oxygen', 'glucose']:
            # Interpolate simulation results to experimental positions
            sim_conc = np.interp(
                exp_data[species]['radial_position'],
                simulation_results[species]['position'],
                simulation_results[species]['concentration']
            )
            
            # Calculate metrics
            error = np.mean(
                np.abs(sim_conc - exp_data[species]['concentration']) /
                exp_data[species]['concentration']
            )
            
            correlation, _ = pearsonr(sim_conc, 
                                    exp_data[species]['concentration'])
            
            metrics[species] = {
                'error': error,
                'correlation': correlation,
                'simulation': sim_conc
            }
            
        self.validation_metrics['biochemical'] = metrics
        return metrics
        
    def validate_cellular_response(self, simulation_results):
        """Validate cellular behavior"""
        self.simulation_results = simulation_results
        exp_data = self.experimental_data['cellular']
        
        # Interpolate simulation viability to experimental oxygen levels
        sim_viability = np.interp(
            exp_data['viability']['oxygen_level'],
            simulation_results['oxygen'],
            simulation_results['viability']
        )
        
        # Calculate metrics
        error = np.mean(
            np.abs(sim_viability - exp_data['viability']['viability']) /
            exp_data['viability']['viability']
        )
        
        correlation, _ = pearsonr(sim_viability,
                                exp_data['viability']['viability'])
        
        self.validation_metrics['cellular'] = {
            'viability_error': error,
            'viability_correlation': correlation,
            'simulation': sim_viability
        }
        
        return self.validation_metrics['cellular']
        
    def run_full_validation(self, simulation_results):
        """Run all validation tests"""
        mechanical_metrics = self.validate_mechanical_response(
            simulation_results['mechanical']
        )
        
        biochemical_metrics = self.validate_biochemical_distribution(
            simulation_results['biochemical']
        )
        
        cellular_metrics = self.validate_cellular_response(
            simulation_results['cellular']
        )
        
        return {
            'mechanical': mechanical_metrics,
            'biochemical': biochemical_metrics,
            'cellular': cellular_metrics
        }
        
    def plot_validation_results(self, save_dir=None):
        """Plot validation comparisons"""
        if self.simulation_results is None:
            raise ValueError("No simulation results available. Run validation first.")
            
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Mechanical response
        ax = axes[0, 0]
        mech_data = self.experimental_data['mechanical']['compression']
        ax.plot(mech_data['strain'], mech_data['stress'], 'o-', 
                label='Experimental', color='blue')
        ax.plot(self.simulation_results['mechanical']['strain'],
                self.simulation_results['mechanical']['stress'], 's--',
                label='Simulation', color='red')
        ax.set_xlabel('Strain')
        ax.set_ylabel('Stress (MPa)')
        ax.set_title('Mechanical Response')
        ax.grid(True)
        ax.legend()
        
        # Oxygen distribution
        ax = axes[0, 1]
        oxygen_data = self.experimental_data['biochemical']['oxygen']
        ax.plot(oxygen_data['radial_position'], 
                oxygen_data['concentration'], 'o-',
                label='Experimental', color='blue')
        if 'simulation' in self.validation_metrics['biochemical']['oxygen']:
            ax.plot(oxygen_data['radial_position'],
                   self.validation_metrics['biochemical']['oxygen']['simulation'],
                   's--', label='Simulation', color='red')
        ax.set_xlabel('Radial Position')
        ax.set_ylabel('Oxygen Concentration')
        ax.set_title('Oxygen Distribution')
        ax.grid(True)
        ax.legend()
        
        # Glucose distribution
        ax = axes[1, 0]
        glucose_data = self.experimental_data['biochemical']['glucose']
        ax.plot(glucose_data['radial_position'],
                glucose_data['concentration'], 'o-',
                label='Experimental', color='blue')
        if 'simulation' in self.validation_metrics['biochemical']['glucose']:
            ax.plot(glucose_data['radial_position'],
                   self.validation_metrics['biochemical']['glucose']['simulation'],
                   's--', label='Simulation', color='red')
        ax.set_xlabel('Radial Position')
        ax.set_ylabel('Glucose Concentration')
        ax.set_title('Glucose Distribution')
        ax.grid(True)
        ax.legend()
        
        # Cell viability
        ax = axes[1, 1]
        viability_data = self.experimental_data['cellular']['viability']
        ax.plot(viability_data['oxygen_level'],
                viability_data['viability'], 'o-',
                label='Experimental', color='blue')
        if 'simulation' in self.validation_metrics['cellular']:
            ax.plot(viability_data['oxygen_level'],
                   self.validation_metrics['cellular']['simulation'],
                   's--', label='Simulation', color='red')
        ax.set_xlabel('Oxygen Level')
        ax.set_ylabel('Cell Viability')
        ax.set_title('Cell Viability vs Oxygen')
        ax.grid(True)
        ax.legend()
        
        plt.tight_layout()
        
        if save_dir:
            save_path = Path(save_dir) / 'validation_plots.png'
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
        
    def generate_validation_report(self, save_dir=None):
        """Generate a comprehensive validation report"""
        if not self.validation_metrics:
            raise ValueError("No validation metrics available. Run validation first.")
            
        report = ["# IVD Model Validation Report\n"]
        
        # Mechanical validation results
        report.append("## Mechanical Validation")
        mech = self.validation_metrics['mechanical']
        report.append(f"- Stress Error: {mech['stress_error']:.2%}")
        report.append(f"- Stress Correlation: {mech['stress_correlation']:.3f}")
        report.append(f"- Displacement Error: {mech['displacement_error']:.2%}\n")
        
        # Biochemical validation results
        report.append("## Biochemical Validation")
        for species in ['oxygen', 'glucose']:
            bio = self.validation_metrics['biochemical'][species]
            report.append(f"\n### {species.capitalize()}")
            report.append(f"- Concentration Error: {bio['error']:.2%}")
            report.append(f"- Correlation: {bio['correlation']:.3f}")
            
        # Cellular validation results
        report.append("\n## Cellular Validation")
        cell = self.validation_metrics['cellular']
        report.append(f"- Viability Error: {cell['viability_error']:.2%}")
        report.append(f"- Viability Correlation: {cell['viability_correlation']:.3f}")
        
        # Write report
        if save_dir:
            report_path = Path(save_dir) / 'validation_report.md'
            with open(report_path, 'w') as f:
                f.write('\n'.join(report))
                
        return '\n'.join(report)
