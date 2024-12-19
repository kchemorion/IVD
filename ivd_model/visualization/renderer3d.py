import plotly.graph_objects as go
import numpy as np
from pathlib import Path
import json

class Renderer3D:
    """Handles 3D visualization of the IVD model"""
    
    def __init__(self, save_dir):
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(exist_ok=True)
        
        # Color maps for different fields
        self.colormaps = {
            'stress': 'Viridis',
            'strain': 'Viridis',
            'pressure': 'RdBu',
            'oxygen': 'YlOrRd',
            'glucose': 'YlOrRd',
            'lactate': 'YlOrRd',
            'ph': 'RdBu',
            'cell_density': 'YlOrRd'
        }
        
        # Units and ranges for different fields
        self.field_info = {
            'stress': {'units': 'MPa', 'range': [-2, 2]},
            'strain': {'units': '-', 'range': [-0.3, 0.3]},
            'pressure': {'units': 'MPa', 'range': [-1, 1]},
            'oxygen': {'units': 'mmHg', 'range': [0, 100]},
            'glucose': {'units': 'mM', 'range': [0, 5]},
            'lactate': {'units': 'mM', 'range': [0, 10]},
            'ph': {'units': '-', 'range': [6.8, 7.4]},
            'cell_density': {'units': 'cells/mm³', 'range': [0, 1e6]}
        }
    
    def create_figure(self, mesh_data, field_name, time=0):
        """Create a 3D figure with field data"""
        vertices = mesh_data['vertices']
        faces = mesh_data['faces']
        values = mesh_data['vertex_values'][field_name]
        
        # Create the mesh3d trace
        fig = go.Figure(data=[
            go.Mesh3d(
                x=vertices[:, 0],
                y=vertices[:, 1],
                z=vertices[:, 2],
                i=faces[:, 0],
                j=faces[:, 1],
                k=faces[:, 2],
                intensity=values,
                colorscale=self.colormaps[field_name],
                colorbar=dict(
                    title=f"{field_name} ({self.field_info[field_name]['units']})"
                ),
                cmin=self.field_info[field_name]['range'][0],
                cmax=self.field_info[field_name]['range'][1]
            )
        ])
        
        # Update layout
        fig.update_layout(
            title=f"{field_name.replace('_', ' ').title()} at t = {time:.1f} days",
            scene=dict(
                aspectmode='data',
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.5),
                    center=dict(x=0, y=0, z=0)
                )
            ),
            width=800,
            height=800
        )
        
        return fig
    
    def save_visualization(self, fig, filename):
        """Save visualization as HTML and capture frame for animation"""
        # Save interactive HTML
        html_path = self.save_dir / f"{filename}.html"
        fig.write_html(str(html_path))
        
        # Save frame data for animation
        frame_data = {
            'data': fig.data[0],
            'layout': fig.layout
        }
        json_path = self.save_dir / f"{filename}_frame.json"
        with open(json_path, 'w') as f:
            json.dump(frame_data, f)
        
    def create_animation(self, frames_pattern, output_name):
        """Create animation from saved frames"""
        # Load all frame data
        frame_files = sorted(self.save_dir.glob(frames_pattern))
        frames = []
        
        for frame_file in frame_files:
            with open(frame_file, 'r') as f:
                frame_data = json.load(f)
                frames.append(frame_data)
        
        if frames:
            # Create figure with animation
            fig = go.Figure(
                data=frames[0]['data'],
                layout=frames[0]['layout']
            )
            
            # Add animation frames
            fig.frames = [go.Frame(data=frame['data']) for frame in frames]
            
            # Add animation buttons
            fig.update_layout(
                updatemenus=[{
                    'buttons': [
                        {
                            'args': [None, {'frame': {'duration': 500, 'redraw': True},
                                          'fromcurrent': True}],
                            'label': 'Play',
                            'method': 'animate'
                        },
                        {
                            'args': [[None], {'frame': {'duration': 0, 'redraw': True},
                                            'mode': 'immediate'}],
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
                }]
            )
            
            # Save animation
            animation_path = self.save_dir / f"{output_name}.html"
            fig.write_html(str(animation_path))
            
            return True
        return False
