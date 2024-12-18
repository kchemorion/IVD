# visualization/plotters.py

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.interpolate import interp1d

class ResultsPlotter:
    """Plot simulation results"""
    
    def __init__(self, save_dir='results'):
        self.save_dir = save_dir
        self.setup_style()
        
    def setup_style(self):
        """Set up plotting style"""
        # Use seaborn defaults instead of style file
        sns.set_theme()
        sns.set_palette("husl")
        plt.rcParams['figure.figsize'] = [10, 6]
        plt.rcParams['figure.dpi'] = 100
        
    def plot_time_series(self, time, data, title, ylabel, filename,
                        uncertainty=None):
        """Plot time series with optional uncertainty bounds"""
        plt.figure()
        
        if uncertainty is not None:
            plt.fill_between(time, 
                           uncertainty['lower'], 
                           uncertainty['upper'],
                           alpha=0.3)
        
        plt.plot(time, data, '-')
        plt.title(title)
        plt.xlabel('Time (days)')
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.savefig(f"{self.save_dir}/{filename}")
        plt.close()
        
    def plot_spatial_distribution(self, x, y, z, title, filename):
        """Plot 2D spatial distribution"""
        plt.figure()
        
        # Ensure z has the same length as x and y
        if len(z) != len(x):
            # Interpolate z to match mesh points if needed
            from scipy.interpolate import griddata
            xi = np.linspace(x.min(), x.max(), len(x))
            yi = np.linspace(y.min(), y.max(), len(y))
            zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
            z = zi.flatten()
        
        # Create triangulation for plotting
        plt.tripcolor(x, y, z, shading='gouraud')
        plt.colorbar(label='Value')
        plt.title(title)
        plt.xlabel('R (mm)')
        plt.ylabel('Z (mm)')
        plt.axis('equal')
        plt.savefig(f"{self.save_dir}/{filename}")
        plt.close()
        
    def plot_sensitivity_analysis(self, sensitivity_results, filename):
        """Plot sensitivity analysis results"""
        plt.figure(figsize=(12, 6))
        
        # Plot first-order indices
        plt.subplot(121)
        plt.bar(range(len(sensitivity_results['S1'])), 
                sensitivity_results['S1'])
        plt.xticks(range(len(sensitivity_results['names'])),
                  sensitivity_results['names'], rotation=45)
        plt.title('First-order Sensitivity Indices')
        
        # Plot total indices
        plt.subplot(122)
        plt.bar(range(len(sensitivity_results['ST'])), 
                sensitivity_results['ST'])
        plt.xticks(range(len(sensitivity_results['names'])),
                  sensitivity_results['names'], rotation=45)
        plt.title('Total Sensitivity Indices')
        
        plt.tight_layout()
        plt.savefig(f"{self.save_dir}/{filename}")
        plt.close()
        
    def plot_validation_comparison(self, predicted, experimental, title, filename):
        """Plot model validation comparison with proper interpolation"""
        plt.figure()
        
        # Debug prints
        print(f"Debugging {title}:")
        print(f"Predicted shape: {np.array(predicted).shape}")
        print(f"Experimental times shape: {np.array(experimental['time_points']).shape}")
        print(f"Experimental values shape: {np.array(experimental['values']).shape}")
        
        # Convert lists to numpy arrays and ensure correct shapes
        exp_times = np.array(experimental['time_points'])
        exp_values = np.array(experimental['values'])
        exp_std = np.array(experimental['std_dev'])
        pred_values = np.array(predicted)
        
        # Make sure pred_times matches pred_values length
        pred_times = np.linspace(0, max(exp_times), len(pred_values))
        
        print(f"pred_times shape: {pred_times.shape}")
        print(f"pred_values shape: {pred_values.shape}")
        
        try:
            # Plot experimental data with error bars
            plt.errorbar(exp_times, 
                        exp_values,
                        yerr=exp_std,
                        fmt='o', 
                        label='Experimental',
                        capsize=5,
                        markersize=8)
            
            # Plot model predictions directly
            plt.plot(pred_times, pred_values, '-', label='Model')
            
            plt.title(title)
            plt.xlabel('Time (days)')
            plt.ylabel(title.split()[0])
            plt.legend()
            plt.grid(True)
            
            # Add confidence interval
            plt.fill_between(exp_times,
                            exp_values - exp_std,
                            exp_values + exp_std,
                            alpha=0.2,
                            label='Experimental ±σ')
            
            plt.tight_layout()
            plt.savefig(f"{self.save_dir}/{filename}")
            plt.close()
            
        except Exception as e:
            print(f"Error in plotting: {str(e)}")