# validation/validation_metrics.py

import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
from scipy.interpolate import interp1d

class ValidationMetrics:
    """Calculate validation metrics for model outputs"""
    
    def calculate_rmse(self, predicted, actual, time_points=None, actual_times=None):
        """Calculate RMSE with safety checks and proper interpolation"""
        try:
            if time_points is not None and actual_times is not None:
                # Create interpolation function for predictions to match experimental time points
                f_pred = interp1d(time_points, predicted, kind='linear', 
                                bounds_error=False, fill_value='extrapolate')
                
                # Interpolate predictions to experimental time points
                pred_interp = f_pred(actual_times)
                
                # Remove any infinite or NaN values
                mask = np.isfinite(pred_interp) & np.isfinite(actual)
                if not np.any(mask):
                    return float('inf')
                
                # Use only valid values
                pred_clean = np.clip(pred_interp[mask], -1e6, 1e6)
                actual_clean = np.clip(actual[mask], -1e6, 1e6)
                
                return np.sqrt(mean_squared_error(actual_clean, pred_clean))
            else:
                return float('inf')
        except:
            return float('inf')
    
    def calculate_r2(self, predicted, actual, time_points=None, actual_times=None):
        """Calculate R2 score with interpolation"""
        try:
            if time_points is not None and actual_times is not None:
                f_pred = interp1d(time_points, predicted, kind='linear', 
                                bounds_error=False, fill_value='extrapolate')
                pred_interp = f_pred(actual_times)
                
                mask = np.isfinite(pred_interp) & np.isfinite(actual)
                if not np.any(mask):
                    return 0.0
                
                pred_clean = np.clip(pred_interp[mask], -1e6, 1e6)
                actual_clean = np.clip(actual[mask], -1e6, 1e6)
                
                return r2_score(actual_clean, pred_clean)
            else:
                return 0.0
        except:
            return 0.0
    
    def calculate_relative_error(self, predicted, actual, time_points=None, actual_times=None):
        """Calculate relative error with interpolation"""
        try:
            if time_points is not None and actual_times is not None:
                f_pred = interp1d(time_points, predicted, kind='linear', 
                                bounds_error=False, fill_value='extrapolate')
                pred_interp = f_pred(actual_times)
                
                mask = np.isfinite(pred_interp) & np.isfinite(actual)
                if not np.any(mask):
                    return float('inf')
                
                pred_clean = np.clip(pred_interp[mask], -1e6, 1e6)
                actual_clean = np.clip(actual[mask], -1e6, 1e6)
                
                # Avoid division by zero
                actual_clean = np.where(actual_clean == 0, 1e-10, actual_clean)
                
                return np.mean(np.abs(pred_clean - actual_clean) / np.abs(actual_clean))
            else:
                return float('inf')
        except:
            return float('inf')
    
    def validate_results(self, predictions, experimental):
        """Validate model predictions against experimental data"""
        results = {}
        
        for metric in predictions:
            if metric in experimental:
                pred = np.array(predictions[metric])
                actual = np.array(experimental[metric]['values'])
                times = np.array(experimental[metric]['time_points'])
                pred_times = np.linspace(0, times[-1], len(pred))
                
                results[metric] = {
                    'rmse': self.calculate_rmse(pred, actual, pred_times, times),
                    'r2': self.calculate_r2(pred, actual, pred_times, times),
                    'relative_error': self.calculate_relative_error(pred, actual, pred_times, times)
                }
                
                # Add threshold check
                rmse_threshold = 0.2  # You can adjust this value
                results[metric]['passes_threshold'] = (
                    results[metric]['rmse'] < rmse_threshold
                )
                
        return results