import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd
from pathlib import Path

class IVDVisualizer:
    """
    Comprehensive visualization tools for IVD model results
    """
    def __init__(self, geometry_model):
        self.geometry = geometry_model
        self.setup_colormaps()
        self.output_dir = Path("visualization_output")
        self.output_dir.mkdir(exist_ok=True)
        
    def setup_colormaps(self):
        """Setup colormaps for different field variables"""
        self.colormaps = {
            'stress': 'RdYlBu_r',      # Red-Yellow-Blue (reversed)
            'strain': 'RdYlBu_r',
            'pressure': 'RdBu',        # Red-Blue
            'oxygen': 'YlOrRd',        # Yellow-Orange-Red
            'glucose': 'YlOrRd',
            'lactate': 'YlOrRd',
            'ph': 'RdBu',
            'cell_density': 'Viridis',
            'cell_viability': 'RdYlGn', # Red-Yellow-Green
            'collagen': 'Magma',
            'proteoglycan': 'Magma',
            'inflammatory': 'Reds'
        }
        
        # Value ranges for each variable
        self.value_ranges = {
            'stress': [-2.0, 2.0],      # MPa
            'strain': [-0.3, 0.3],      # dimensionless
            'pressure': [-0.5, 0.5],    # MPa
            'oxygen': [0, 7.0],         # kPa
            'glucose': [0, 5.0],        # mM
            'lactate': [0, 10.0],       # mM
            'ph': [6.8, 7.4],           # pH units
            'cell_density': [0, 8e6],   # cells/ml
            'cell_viability': [0, 1],   # fraction
            'collagen': [0, 0.2],       # g/ml
            'proteoglycan': [0, 0.1],   # g/ml
            'inflammatory': [0, 1]       # normalized
        }
        
    def visualize_3d_field(self, field_name, field_data, time=0):
        """
        Create 3D visualization of a field variable on the IVD geometry
        """
        # Get mesh data
        vertices = self.geometry.mesh.points
        faces = self.geometry.mesh.cells[0].data
        
        # Create figure
        fig = go.Figure()
        
        # Add mesh surface with field data
        fig.add_trace(go.Mesh3d(
            x=vertices[:, 0],
            y=vertices[:, 1],
            z=vertices[:, 2],
            i=faces[:, 0],
            j=faces[:, 1],
            k=faces[:, 2],
            intensity=field_data,
            colorscale=self.colormaps[field_name],
            cmin=self.value_ranges[field_name][0],
            cmax=self.value_ranges[field_name][1],
            colorbar=dict(
                title=field_name,
                tickformat='.2f'
            ),
            lighting=dict(
                ambient=0.5,
                diffuse=0.8,
                specular=0.1,
                roughness=0.1
            ),
            lightposition=dict(
                x=100,
                y=200,
                z=150
            )
        ))
        
        # Update layout
        fig.update_layout(
            title=f"{field_name.replace('_', ' ').title()} at t = {time:.1f} days",
            scene=dict(
                aspectmode='data',
                camera=dict(
                    up=dict(x=0, y=0, z=1),
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=1.5, y=1.5, z=1.5)
                ),
                xaxis_title='X (mm)',
                yaxis_title='Y (mm)',
                zaxis_title='Z (mm)'
            ),
            width=800,
            height=800
        )
        
        return fig
        
    def create_field_animation(self, field_name, field_data_series, times):
        """
        Create animation of field evolution over time
        """
        # Create figure with initial state
        fig = self.visualize_3d_field(field_name, field_data_series[0], times[0])
        
        # Add frames for animation
        frames = []
        for i, (field_data, time) in enumerate(zip(field_data_series, times)):
            frame = go.Frame(
                data=[go.Mesh3d(
                    intensity=field_data,
                    colorscale=self.colormaps[field_name],
                    cmin=self.value_ranges[field_name][0],
                    cmax=self.value_ranges[field_name][1]
                )],
                name=f'frame{i}'
            )
            frames.append(frame)
            
        fig.frames = frames
        
        # Add animation controls
        fig.update_layout(
            updatemenus=[{
                'buttons': [
                    {
                        'args': [None, {'frame': {'duration': 500, 'redraw': True},
                                      'fromcurrent': True,
                                      'transition': {'duration': 0}}],
                        'label': 'Play',
                        'method': 'animate'
                    },
                    {
                        'args': [[None], {'frame': {'duration': 0, 'redraw': True},
                                        'mode': 'immediate',
                                        'transition': {'duration': 0}}],
                        'label': 'Pause',
                        'method': 'animate'
                    }
                ],
                'direction': 'left',
                'pad': {'r': 10, 't': 87},
                'showactive': False,
                'type': 'buttons',
                'x': 0.1,
                'xanchor': 'right',
                'y': 0,
                'yanchor': 'top'
            }],
            sliders=[{
                'currentvalue': {
                    'font': {'size': 20},
                    'prefix': 'Time: ',
                    'suffix': ' days',
                    'visible': True,
                    'xanchor': 'right'
                },
                'pad': {'b': 10, 't': 50},
                'len': 0.9,
                'x': 0.1,
                'y': 0,
                'steps': [
                    {
                        'args': [[f'frame{i}'],
                                {'frame': {'duration': 500, 'redraw': True},
                                 'mode': 'immediate',
                                 'transition': {'duration': 0}}],
                        'label': f'{time:.1f}',
                        'method': 'animate'
                    }
                    for i, time in enumerate(times)
                ]
            }]
        )
        
        return fig
        
    def create_cross_section_view(self, field_name, field_data, plane='sagittal', position=0.5):
        """
        Create cross-sectional view of the field
        
        Parameters:
        -----------
        plane : str
            'sagittal', 'coronal', or 'axial'
        position : float
            Position of the cutting plane (0-1)
        """
        vertices = self.geometry.mesh.points
        
        # Get plane coordinates
        if plane == 'sagittal':
            x_range = vertices[:, 0].max() - vertices[:, 0].min()
            cut_pos = vertices[:, 0].min() + position * x_range
            mask = np.abs(vertices[:, 0] - cut_pos) < 0.5  # mm tolerance
            plot_x, plot_y = vertices[mask, 1], vertices[mask, 2]
        elif plane == 'coronal':
            y_range = vertices[:, 1].max() - vertices[:, 1].min()
            cut_pos = vertices[:, 1].min() + position * y_range
            mask = np.abs(vertices[:, 1] - cut_pos) < 0.5
            plot_x, plot_y = vertices[mask, 0], vertices[mask, 2]
        else:  # axial
            z_range = vertices[:, 2].max() - vertices[:, 2].min()
            cut_pos = vertices[:, 2].min() + position * z_range
            mask = np.abs(vertices[:, 2] - cut_pos) < 0.5
            plot_x, plot_y = vertices[mask, 0], vertices[mask, 1]
            
        # Create scatter plot
        fig = go.Figure(data=go.Scatter(
            x=plot_x,
            y=plot_y,
            mode='markers',
            marker=dict(
                size=5,
                color=field_data[mask],
                colorscale=self.colormaps[field_name],
                cmin=self.value_ranges[field_name][0],
                cmax=self.value_ranges[field_name][1],
                colorbar=dict(title=field_name)
            )
        ))
        
        # Update layout
        fig.update_layout(
            title=f"{field_name.replace('_', ' ').title()} - {plane.title()} View",
            xaxis_title='Position (mm)',
            yaxis_title='Position (mm)',
            width=600,
            height=600
        )
        
        return fig
        
    def create_time_series_plot(self, time_series_data):
        """
        Create time series plots for multiple variables
        """
        # Create subplots
        n_vars = len(time_series_data)
        fig = make_subplots(rows=n_vars, cols=1,
                           subplot_titles=[key.replace('_', ' ').title() 
                                         for key in time_series_data.keys()])
        
        # Add traces for each variable
        for i, (var_name, data) in enumerate(time_series_data.items(), 1):
            fig.add_trace(
                go.Scatter(x=data['time'], y=data['value'],
                          name=var_name.replace('_', ' ').title(),
                          line=dict(width=2)),
                row=i, col=1
            )
            
            # Add experimental data if available
            if 'experimental' in data:
                fig.add_trace(
                    go.Scatter(x=data['experimental']['time'],
                              y=data['experimental']['value'],
                              name=f'{var_name} (Experimental)',
                              mode='markers',
                              marker=dict(size=8)),
                    row=i, col=1
                )
                
        # Update layout
        fig.update_layout(
            height=300 * n_vars,
            width=800,
            showlegend=True,
            title_text="Time Evolution of IVD Properties"
        )
        
        return fig

    def create_interactive_dashboard(self, simulation_results):
        """
        Create an interactive dashboard with multiple visualizations
        """
        # Extract time points
        times = simulation_results['time']
        
        # Create subplot figure
        fig = make_subplots(
            rows=2, cols=2,
            specs=[[{'type': 'mesh3d'}, {'type': 'scatter'}],
                   [{'type': 'scatter3d'}, {'type': 'scatter'}]],
            subplot_titles=('3D Visualization', 'Time Series',
                          'Cross Sections', 'Property Correlations')
        )
        
        # Add 3D visualization
        vertices = self.geometry.mesh.points
        faces = self.geometry.mesh.cells[0].data
        current_field = simulation_results['fields']['stress'][-1]
        
        fig.add_trace(
            go.Mesh3d(
                x=vertices[:, 0],
                y=vertices[:, 1],
                z=vertices[:, 2],
                i=faces[:, 0],
                j=faces[:, 1],
                k=faces[:, 2],
                intensity=current_field,
                colorscale=self.colormaps['stress'],
                colorbar=dict(title='Stress (MPa)')
            ),
            row=1, col=1
        )
        
        # Add time series
        for field_name, field_data in simulation_results['fields'].items():
            fig.add_trace(
                go.Scatter(
                    x=times,
                    y=np.mean(field_data, axis=1),
                    name=field_name.replace('_', ' ').title()
                ),
                row=1, col=2
            )
            
        # Add cross sections
        z_levels = np.linspace(vertices[:, 2].min(),
                             vertices[:, 2].max(), 5)
        for z in z_levels:
            mask = np.abs(vertices[:, 2] - z) < 0.5
            fig.add_trace(
                go.Scatter3d(
                    x=vertices[mask, 0],
                    y=vertices[mask, 1],
                    z=vertices[mask, 2],
                    mode='markers',
                    marker=dict(
                        size=3,
                        color=current_field[mask],
                        colorscale=self.colormaps['stress']
                    ),
                    showlegend=False
                ),
                row=2, col=1
            )
            
        # Add correlations
        field1 = simulation_results['fields']['stress'][-1]
        field2 = simulation_results['fields']['strain'][-1]
        fig.add_trace(
            go.Scatter(
                x=field1,
                y=field2,
                mode='markers',
                marker=dict(
                    color=simulation_results['fields']['oxygen'][-1],
                    colorscale=self.colormaps['oxygen'],
                    showscale=True,
                    colorbar=dict(title='Oxygen (kPa)')
                ),
                name='Stress-Strain'
            ),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            height=1000,
            width=1200,
            title_text="IVD Model Analysis Dashboard",
            showlegend=True
        )
        
        return fig

    def save_visualization(self, fig, filename):
        """Save visualization to file"""
        filepath = self.output_dir / filename
        fig.write_html(str(filepath))