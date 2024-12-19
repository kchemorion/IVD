    def update_growth_factors(self, treatment, dt, effects):
        """Update growth factor diffusion and degradation"""
        state = treatment['state']
        
        for factor, concentration in state.items():
            # Get diffusion coefficient
            D = self.diffusion_coef['growth_factors'][factor]
            
            # Solve diffusion equation
            new_conc = self.solve_diffusion(concentration, D, dt)
            
            # Apply degradation
            degradation_rate = self.params.GROWTH_FACTOR_DEGRADATION[factor]
            new_conc *= (1 - degradation_rate * dt)
            
            # Update state and effects
            state[factor] = new_conc
            effects['growth_factors'][factor] += new_conc
            
    def update_gene_therapy(self, treatment, dt, effects):
        """Update gene therapy treatment"""
        state = treatment['state']
        
        # Update transfection efficiency
        decay_rate = self.params.GENE_THERAPY_DECAY
        state['transfection'] *= (1 - decay_rate * dt)
        
        # Update gene expression
        expression_rate = self.params.GENE_EXPRESSION_RATE
        target_expression = state['transfection']
        current_expression = state['expression']
        
        d_expression = expression_rate * (target_expression - current_expression) * dt
        state['expression'] += d_expression
        
        # Add to effects
        effects['gene_expression'] += state['expression']
        
    def update_biomaterial(self, treatment, dt, mechanical_state, effects):
        """Update biomaterial degradation and mechanics"""
        state = treatment['state']
        
        # Update crosslinking
        crosslink_rate = self.params.BIOMATERIAL_CROSSLINK_RATE
        max_crosslinks = state['concentration']  # Maximum possible crosslinks
        current_crosslinks = state['crosslinking']
        
        d_crosslinks = crosslink_rate * (max_crosslinks - current_crosslinks) * dt
        state['crosslinking'] += d_crosslinks
        
        # Update degradation
        strain_energy = np.sum(mechanical_state['strain']**2)
        degradation_rate = (self.params.BIOMATERIAL_BASE_DEGRADATION * 
                          (1 + self.params.BIOMATERIAL_MECH_DEGRADATION * strain_energy))
        
        state['concentration'] *= (1 - degradation_rate * dt)
        state['crosslinking'] *= (1 - degradation_rate * dt)
        
        # Add to effects
        effects['biomaterial'] += state['concentration']
        
    def solve_diffusion(self, concentration, diffusion_coef, dt):
        """Solve diffusion equation for one time step"""
        # Get mesh data
        nodes = self.nodes
        elements = self.elements
        
        # Create FEM matrices (simplified)
        n_nodes = len(nodes)
        M = np.eye(n_nodes)  # Mass matrix
        K = np.zeros((n_nodes, n_nodes))  # Stiffness matrix
        
        # Assemble stiffness matrix
        for elem in elements:
            # Calculate element matrices
            Ke = self.element_diffusion_matrix(nodes[elem], diffusion_coef)
            
            # Add to global matrix
            for i, ii in enumerate(elem):
                for j, jj in enumerate(elem):
                    K[ii, jj] += Ke[i, j]
        
        # Solve system
        A = M/dt + K
        b = M/dt @ concentration
        
        # Add boundary conditions here if needed
        
        # Solve
        new_concentration = np.linalg.solve(A, b)
        return new_concentration
        
    def element_diffusion_matrix(self, nodes, D):
        """Calculate element diffusion matrix"""
        # Calculate element volume
        v1 = nodes[1] - nodes[0]
        v2 = nodes[2] - nodes[0]
        v3 = nodes[3] - nodes[0]
        volume = abs(np.dot(v1, np.cross(v2, v3))) / 6
        
        # Shape function gradients (simplified)
        B = np.array([
            [-1, -1, -1],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])
        
        # Element matrix
        Ke = D * volume * (B @ B.T)
        return Ke
        
    def get_treatment_state(self):
        """Get current state of all treatments"""
        state = {
            'active_treatments': len(self.active_treatments),
            'treatments': {}
        }
        
        for i, treatment in enumerate(self.active_treatments):
            state['treatments'][f'treatment_{i}'] = {
                'type': treatment['type'],
                'time': treatment['time'],
                'state': {k: v.copy() if isinstance(v, np.ndarray) else v 
                         for k, v in treatment['state'].items()}
            }
            
        return state