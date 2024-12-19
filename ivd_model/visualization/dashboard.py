from dash import Dash, html, dcc, Input, Output, State, callback
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from pathlib import Path
from .visualizer import IVDVisualizer

class IVDDashboard:
    """
    Interactive dashboard for IVD model visualization and analysis
    """
    def __init__(self, geometry_model, simulation_runner):
        self.geometry = geometry_model
        self.simulation = simulation_runner
        self.visualizer = IVDVisualizer(geometry_model)
        self.app = Dash(__name__)
        self.setup_layout()
        self.setup_callbacks()
        
    def setup_layout(self):
        """Setup dashboard layout"""
        self.app.layout = html.Div([
            # Header
            html.H1("Intervertebral Disc Model Dashboard",
                   style={'textAlign': 'center'}),
            
            # Main content
            html.Div([
                # Left panel - Controls
                html.Div([
                    html.H3("Visualization Controls"),
                    
                    # Field selection
                    html.Label("Field Variable:"),
                    dcc.Dropdown(
                        id='field-selector',
                        options=[
                            {'label': 'Stress', 'value': 'stress'},
                            {'label': 'Strain', 'value': 'strain'},
                            {'label': 'Pressure', 'value': 'pressure'},
                            {'label': 'Oxygen', 'value': 'oxygen'},
                            {'label': 'Glucose', 'value': 'glucose'},
                            {'label': 'Cell Density', 'value': 'cell_density'},
                            {'label': 'Cell Viability', 'value': 'cell_viability'}
                        ],
                        value='stress'
                    ),
                    
                    # View controls
                    html.Label("View Type:"),
                    dcc.RadioItems(
                        id='view-selector',
                        options=[
                            {'label': '3D View', 'value': '3d'},
                            {'label': 'Sagittal Section', 'value': 'sagittal'},
                            {'label': 'Coronal Section', 'value': 'coronal'},
                            {'label': 'Axial Section', 'value': 'axial'}
                        ],
                        value='3d'
                    ),
                    
                    # Time controls
                    html.Label("Time Point:"),
                    dcc.Slider(
                        id='time-slider',
                        min=0,
                        max=100,
                        step=1,
                        value=0,
                        marks={i: f'{i}d' for i in range(0, 101, 10)}
                    ),
                    
                    # Analysis controls
                    html.H3("Analysis Controls"),
                    
                    # Loading conditions
                    html.Label("Loading Condition:"),
                    dcc.Dropdown(
                        id='loading-selector',
                        options=[
                            {'label': 'Daily Activity', 'value': 'daily'},
                            {'label': 'Exercise', 'value': 'exercise'},
                            {'label': 'Static Load', 'value': 'static'}
                        ],
                        value='daily'
                    ),
                    
                    # Treatment options
                    html.Label("Treatment:"),
                    dcc.Dropdown(
                        id='treatment-selector',
                        options=[
                            {'label': 'None', 'value': 'none'},
                            {'label': 'Cell Therapy', 'value': 'cell'},
                            {'label': 'Growth Factors', 'value': 'growth_factors'},
                            {'label': 'Biomaterial', 'value': 'biomaterial'}
                        ],
                        value='none'
                    ),
                    
                    # Run simulation button
                    html.Button('Run Simulation', id='run-simulation',
                              style={'marginTop': '20px'})
                    
                ], style={'width': '25%', 'float': 'left', 'padding': '20px'}),
                
                # Right panel - Visualizations
                html.Div([
                    # Main visualization
                    dcc.Graph(id='main-visualization',
                             style={'height': '600px'}),
                    
                    # Time series plots
                    html.Div([
                        dcc.Graph(id='time-series',
                                 style={'height': '300px'})
                    ]),
                    
                    # Additional metrics
                    html.Div([
                        html.H3("Key Metrics"),
                        html.Div(id='metrics-display',
                                style={'padding': '10px'})
                    ])
                ], style={'width': '75%', 'float': 'right', 'padding': '20px'})
            ])
        ])
        
    def setup_callbacks(self):
        """Setup dashboard callbacks"""
        @self.app.callback(
            [Output('main-visualization', 'figure'),
             Output('time-series', 'figure'),
             Output('metrics-display', 'children')],
            [Input('field-selector', 'value'),
             Input('view-selector', 'value'),
             Input('time-slider', 'value'),
             Input('run-simulation', 'n_clicks')],
            [State('loading-selector', 'value'),
             State('treatment-selector', 'value')]
        )
        def update_visualization(field, view, time, n_clicks, loading, treatment):
            # Run simulation if requested
            if n_clicks:
                results = self.run_simulation(loading, treatment)
            else:
                results = self.get_current_results()
                
            # Create main visualization
            if view == '3d':
                main_fig = self.visualizer.visualize_3d_field(
                    field,
                    results['fields'][field][time],
                    time
                )
            else:
                main_fig = self.visualizer.create_cross_section_view(
                    field,
                    results['fields'][field][time],
                    view
                )
                
            # Create time series
            time_series = self.visualizer.create_time_series_plot({
                field: {
                    'time': results['time'],
                    'value': np.mean(results['fields'][field], axis=1)
                } for field in results['fields']
            })
            
            # Calculate metrics
            metrics = self.calculate_metrics(results, time)
            metrics_display = self.create_metrics_display(metrics)
            
            return main_fig, time_series, metrics_display
            
    def run_simulation(self, loading_type, treatment_type):
        """Run simulation with specified conditions"""
        # Setup simulation parameters
        if loading_type == 'daily':
            loading = self.simulation.generate_daily_loads()
        elif loading_type == 'exercise':
            loading = self.simulation.generate_exercise_loads()
        else:
            loading = lambda t: 0.5  # Static load
            
        # Setup treatment if selected
        if treatment_type != 'none':
            treatment_params = self.get_treatment_params(treatment_type)
            self.simulation.treatment.add_treatment(
                treatment_type,
                treatment_params,
                time=0
            )
            
        # Run simulation
        results = self.simulation.run_daily_activity()
        return results
        
    def calculate_metrics(self, results, time):
        """Calculate key metrics from results"""
        metrics = {
            'height_change': (results['mechanical']['height'][time] /
                            results['mechanical']['height'][0] - 1) * 100,
            'mean_cell_viability': np.mean(results['cellular']['viability'][time]),
            'min_oxygen': np.min(results['biochemical']['oxygen'][time]),
            'mean_stress': np.mean(np.abs(results['mechanical']['stress'][time]))
        }
        return metrics
        
    def create_metrics_display(self, metrics):
        """Create HTML display of metrics"""
        return html.Div([
            html.P(f"Height Change: {metrics['height_change']:.1f}%"),
            html.P(f"Cell Viability: {metrics['mean_cell_viability']:.1%}"),
            html.P(f"Minimum Oxygen: {metrics['min_oxygen']:.1f} kPa"),
            html.P(f"Mean Stress: {metrics['mean_stress']:.2f} MPa")
        ])
        
    def get_treatment_params(self, treatment_type):
        """Get default parameters for different treatments"""
        if treatment_type == 'cell':
            return {
                'cell_count': 1e6,
                'injection_sites': [[0, 0, 0]],
                'spread': 0.005
            }
        elif treatment_type == 'growth_factors':
            return {
                'factors': {
                    'tgf_beta': 10.0,
                    'bmp2': 100.0
                },
                'injection_sites': [[0, 0, 0]],
                'spread': 0.005
            }
        elif treatment_type == 'biomaterial':
            return {
                'material_type': 'hydrogel',
                'volume': 0.5,
                'injection_sites': [[0, 0, 0]],
                'spread': 0.005
            }
            
    def run_server(self, debug=True):
        """Run the dashboard server"""
        self.app.run_server(debug=debug)