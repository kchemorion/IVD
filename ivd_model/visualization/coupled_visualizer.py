import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class CoupledVisualizer:
    """
    Visualization tools for coupled IVD model behavior.
    Provides interactive 3D visualization and time series plots.
    """
    
    def __init__(self, coupling_model):
        self.coupling = coupling_model
        self.fig = None
        self.time_series = {
            'time': [],
            'mechanical': {},
            'biochemical': {},
            'cellular': {}
        }
        
    def create_visualization(self):
        """Create main visualization dashboard"""
        self.fig = make_subplots(
            rows=2, cols=2,
            specs=[[{'type': 'surface', 'rowspan': 2}, {'type': 'scatter'}],
                  [None, {'type': 'scatter'}]],
            subplot_titles=('3D Visualization', 'Mechanical Response',
                          'Biochemical Distribution')
        )
        
        # Add 3D geometry
        self.update_3d_visualization()
        
        # Add time series plots
        self.update_time_series()
        
        # Update layout
        self.fig.update_layout(
            height=800,
            title_text="Coupled IVD Model Visualization",
            showlegend=True
        )
        
    def update_3d_visualization(self):
        """Update 3D geometry visualization"""
        # Get mesh data
        nodes = self.coupling.spatial.nodes
        elements = self.coupling.spatial.elements
        regions = self.coupling.spatial.element_regions
        
        # Create mesh visualization
        x, y, z = nodes.T
        i, j, k = elements.T
        
        # Color mapping for different regions
        region_colors = {
            'AF': 'rgb(255,127,14)',
            'NP': 'rgb(44,160,44)',
            'CEP': 'rgb(214,39,40)'
        }
        
        colors = [region_colors[region] for region in regions]
        
        # Add mesh to figure
        self.fig.add_trace(
            go.Mesh3d(
                x=x, y=y, z=z,
                i=i, j=j, k=k,
                color='lightgray',
                opacity=0.8,
                name='IVD Geometry'
            ),
            row=1, col=1
        )
        
        # Add scalar field visualization
        field = self.coupling.spatial.state['stress'][:, 0]  # von Mises stress
        self.add_scalar_field(field, 'Von Mises Stress')
        
    def add_scalar_field(self, field, name):
        """Add scalar field visualization"""
        nodes = self.coupling.spatial.nodes
        
        # Create colormap
        self.fig.add_trace(
            go.Scatter3d(
                x=nodes[:, 0],
                y=nodes[:, 1],
                z=nodes[:, 2],
                mode='markers',
                marker=dict(
                    size=5,
                    color=field,
                    colorscale='Viridis',
                    showscale=True
                ),
                name=name
            ),
            row=1, col=1
        )
        
    def update_time_series(self):
        """Update time series plots"""
        time = self.time_series['time']
        
        # Mechanical response
        if 'stress' in self.time_series['mechanical']:
            self.fig.add_trace(
                go.Scatter(
                    x=time,
                    y=self.time_series['mechanical']['stress'],
                    name='Max Stress',
                    line=dict(color='red')
                ),
                row=1, col=2
            )
            
        # Biochemical distribution
        if 'oxygen' in self.time_series['biochemical']:
            self.fig.add_trace(
                go.Scatter(
                    x=time,
                    y=self.time_series['biochemical']['oxygen'],
                    name='Mean Oxygen',
                    line=dict(color='blue')
                ),
                row=2, col=2
            )
            
    def update_state(self, t):
        """Update visualization with current state"""
        # Store time series data
        self.time_series['time'].append(t)
        
        # Mechanical data
        stress = self.coupling.spatial.state['stress']
        self.time_series['mechanical']['stress'] = np.max(stress)
        
        # Biochemical data
        oxygen = self.coupling.biochemical.state['oxygen']
        self.time_series['biochemical']['oxygen'] = np.mean(oxygen)
        
        # Update plots
        self.update_3d_visualization()
        self.update_time_series()
        
    def create_animation(self, times, states):
        """Create animation of model evolution"""
        import plotly.express as px
        
        frames = []
        for t, state in zip(times, states):
            # Update model state
            self.coupling.spatial.state = state
            
            # Create frame
            self.update_3d_visualization()
            frames.append(go.Frame(
                data=self.fig.data,
                name=f't_{t}'
            ))
            
        # Add slider
        self.fig.update_layout(
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
                'type': 'buttons'
            }],
            sliders=[{
                'currentvalue': {'prefix': 'Time: '},
                'steps': [{'args': [[f't_{t}'], {'frame': {'duration': 0, 'redraw': True}}],
                          'label': f'{t:.1f}',
                          'method': 'animate'} for t in times]
            }]
        )
        
    def save_visualization(self, filename):
        """Save visualization to HTML file"""
        self.fig.write_html(filename)
        
    def display(self):
        """Display visualization in notebook or browser"""
        self.fig.show()