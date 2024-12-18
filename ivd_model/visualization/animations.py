# visualization/animations.py

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.animation as animation

class AnimationGenerator:
    """Generate animations of model results"""
    
    def __init__(self, save_dir='results'):
        self.save_dir = save_dir
        
    def create_3d_ivd_animation(self, spatial_model, biochemical_results, times, filename='ivd_3d.gif'):
        """
        Create 3D animation of IVD with biochemical overlays
        
        Parameters:
        -----------
        spatial_model : SpatialModel
            Model containing mesh and deformation data
        biochemical_results : list
            List of biochemical states over time
        times : array-like
            Time points for animation
        filename : str
            Output filename
        """
        # Create figure for 3D plot
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Get mesh data
        nodes = spatial_model.nodes
        elements = spatial_model.elements
        
        # Function to update the animation
        def update(frame):
            ax.clear()
            
            # Get current state
            current_disp = spatial_model.results['displacement'][frame]
            current_o2 = biochemical_results[frame]['oxygen']
            
            # Update node positions with displacement
            deformed_nodes = nodes + current_disp.reshape(-1, 2)
            
            # Create cylindrical mesh by rotating points
            theta = np.linspace(0, 2*np.pi, 36)
            X = np.outer(deformed_nodes[:,0], np.cos(theta))
            Y = np.outer(deformed_nodes[:,0], np.sin(theta))
            Z = np.repeat(deformed_nodes[:,1][:,np.newaxis], len(theta), axis=1)
            
            # Color based on oxygen concentration
            o2_colors = plt.cm.viridis(current_o2/np.max(current_o2))
            
            # Plot surface
            surf = ax.plot_surface(X, Y, Z, facecolors=o2_colors, alpha=0.8)
            
            # Add colorbar
            if frame == 0:
                fig.colorbar(surf, label='Oxygen Concentration')
            
            # Set labels and title
            ax.set_xlabel('X (mm)')
            ax.set_ylabel('Y (mm)')
            ax.set_zlabel('Z (mm)')
            ax.set_title(f'Time: {times[frame]:.1f} days')
            
            # Set consistent view
            ax.view_init(elev=20, azim=frame)  # Rotate view
            ax.set_box_aspect([1,1,1.2])  # Keep aspect ratio
            
        # Create animation
        anim = FuncAnimation(fig, update,
                           frames=len(times),
                           interval=50)
        
        # Save animation
        writer = PillowWriter(fps=30)
        anim.save(f"{self.save_dir}/{filename}", writer=writer)
        plt.close()
        
    def create_cross_section_animation(self, spatial_model, biochemical_results, times, 
                                     filename='cross_section.gif'):
        """Create animation of IVD cross-section with biochemical overlays"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        def update(frame):
            ax1.clear()
            ax2.clear()
            
            # Get current state
            current_disp = spatial_model.results['displacement'][frame]
            current_o2 = biochemical_results[frame]['oxygen']
            current_strain = spatial_model.results['strain'][frame]
            
            # Plot deformed mesh
            deformed_nodes = spatial_model.nodes + current_disp
            
            # Plot oxygen distribution
            im1 = ax1.tricontourf(deformed_nodes[:,0], deformed_nodes[:,1], 
                                current_o2, levels=20, cmap='viridis')
            ax1.set_title('Oxygen Distribution')
            plt.colorbar(im1, ax=ax1)
            
            # Plot strain distribution
            strain_mag = np.sqrt(np.sum(current_strain**2, axis=1))
            im2 = ax2.tricontourf(deformed_nodes[:,0], deformed_nodes[:,1],
                                strain_mag, levels=20, cmap='coolwarm')
            ax2.set_title('Strain Magnitude')
            plt.colorbar(im2, ax=ax2)
            
            # Set labels and title
            for ax in [ax1, ax2]:
                ax.set_xlabel('R (mm)')
                ax.set_ylabel('Z (mm)')
                ax.set_aspect('equal')
            
            fig.suptitle(f'Time: {times[frame]:.1f} days')
            
        # Create animation
        anim = FuncAnimation(fig, update,
                           frames=len(times),
                           interval=50)
        
        # Save animation
        writer = PillowWriter(fps=30)
        anim.save(f"{self.save_dir}/{filename}", writer=writer)
        plt.close()