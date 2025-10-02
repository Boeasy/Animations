""" Visualizing the stress strain curve of gold (or a general metal) up to the fracture while stretching, visualizing the graph and the atomistic behavior"""

from manim import *
import numpy as np

#So far this thing is in the right shape, it just needs to have the 
#fracture point synced up, maybe exaggerate the movement of the lattice as well.
#It could also help to make the graph a bit more representative of the forces, instead of just a piecewise linear function
#As a stop gap we could just make the second half more steep downwards...
# If time allows, add a section to this to explain the moduli (young's, bulk, shear, poisson)


class StressFracture(Scene):
    def construct(self):
        # Parameters
        rows = 6
        cols = 8
        initial_spacing = 0.5
        atom_radius = 0.15
        
        # Create title
        title = Text("Stress-Strain Behavior: Gold Lattice Fracture", font_size=32)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Create split screen divider
        divider = Line(ORIGIN, 4*DOWN, color=WHITE, stroke_width=2)
        divider.shift(0.3*RIGHT)
        
        # Left side: Lattice
        # Right side: Graph
        
        # Create initial gold lattice
        lattice_atoms = VGroup()
        for i in range(rows):
            for j in range(cols):
                atom = Circle(radius=atom_radius, color=GOLD, fill_opacity=1)
                x = (j - cols/2 + 0.5) * initial_spacing
                y = (i - rows/2 + 0.5) * initial_spacing - 0.5
                atom.move_to([x, y, 0])
                lattice_atoms.add(atom)
        
        # Shift lattice to left side
        lattice_atoms.shift(3*LEFT)
        
        # Create bonds between atoms
        bonds = VGroup()
        for i in range(rows):
            for j in range(cols):
                atom_index = i * cols + j
                atom = lattice_atoms[atom_index]
                
                # Horizontal bonds
                if j < cols - 1:
                    neighbor = lattice_atoms[atom_index + 1]
                    bond = Line(atom.get_center(), neighbor.get_center(), 
                               color=GRAY, stroke_width=2)
                    bonds.add(bond)
                
                # Vertical bonds
                if i < rows - 1:
                    neighbor = lattice_atoms[atom_index + cols]
                    bond = Line(atom.get_center(), neighbor.get_center(), 
                               color=GRAY, stroke_width=2)
                    bonds.add(bond)
        
        # Set up axes for stress-strain curve
        axes = Axes(
            x_range=[0, 0.25, 0.05],
            y_range=[0, 1.2, 0.2],
            x_length=4,
            y_length=3.5,
            axis_config={"include_tip": True, "numbers_to_exclude": []},
            tips=True,
        ).shift(2.5*RIGHT + 0.3*DOWN)
        
        # Add axis labels
        x_label = axes.get_x_axis_label(Text("Strain", font_size=24))
        y_label = axes.get_y_axis_label(Text("Stress", font_size=24), edge=LEFT, direction=LEFT)
        
        # Create stress-strain curve - linear up to fracture, then sharp drop
        fracture_strain = 0.15
        max_strain = 0.20  # Point where material breaks
        
        def stress_strain_curve(strain):
            if strain < fracture_strain:  # Linear elastic region
                return 5.0 * strain  # Linear slope
            else:  # Sharp drop after fracture
                peak_stress = 5.0 * fracture_strain
                # Nearly vertical drop
                return peak_stress * (1 - 10 * (strain - fracture_strain))
        
        curve = axes.plot(stress_strain_curve, x_range=[0, max_strain], color=BLUE, stroke_width=4)
        
        # Mark fracture point on the curve
        fracture_point = axes.coords_to_point(fracture_strain, stress_strain_curve(fracture_strain))
        fracture_dot = Dot(fracture_point, color=RED, radius=0.08)
        fracture_label = Text("Fracture", font_size=20, color=RED).next_to(fracture_dot, UP+RIGHT, buff=0.1)
        
        # Show initial setup
        self.play(
            Create(divider),
            Create(bonds),
            Create(lattice_atoms),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Show the graph first with the complete curve
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            run_time=1.5
        )
        self.wait(0.5)
        
        self.play(
            Create(curve),
            Create(fracture_dot),
            Write(fracture_label),
            run_time=2
        )
        self.wait(0.5)
        
        # Create indicator dot that will follow the curve
        moving_dot = Dot(axes.coords_to_point(0, 0), color=YELLOW, radius=0.08)
        self.add(moving_dot)
        
        # Create strain tracker
        strain = ValueTracker(0)
        total_animation_time = 8.0
        
        # Updater for lattice stretching
        def update_lattice(mob):
            current_strain = strain.get_value()
            
            for i in range(rows):
                for j in range(cols):
                    atom_index = i * cols + j
                    atom = mob[atom_index]
                    
                    # Original position
                    x_orig = (j - cols/2 + 0.5) * initial_spacing
                    y_orig = (i - rows/2 + 0.5) * initial_spacing - 0.5
                    
                    # Apply horizontal strain from center
                    distance_from_center = x_orig
                    x_new = x_orig * (1 + current_strain)
                    
                    # Add some vertical contraction (Poisson effect)
                    poisson_ratio = 0.44  # Gold's Poisson ratio
                    y_new = y_orig * (1 - poisson_ratio * current_strain * 0.5)
                    
                    atom.move_to([x_new - 3, y_new, 0])
        
        # Updater for bonds
        def update_bonds(mob):
            bond_index = 0
            current_strain = strain.get_value()
            
            # Calculate opacity based on strain (bonds weaken as we approach fracture)
            if current_strain < fracture_strain:
                bond_opacity = 1.0
            else:
                bond_opacity = max(0, 1.0 - 5*(current_strain - fracture_strain))
            
            for i in range(rows):
                for j in range(cols):
                    atom_index = i * cols + j
                    atom = lattice_atoms[atom_index]
                    
                    # Horizontal bonds
                    if j < cols - 1:
                        neighbor = lattice_atoms[atom_index + 1]
                        bond = mob[bond_index]
                        bond.put_start_and_end_on(atom.get_center(), neighbor.get_center())
                        bond.set_opacity(bond_opacity)
                        bond_index += 1
                    
                    # Vertical bonds
                    if i < rows - 1:
                        neighbor = lattice_atoms[atom_index + cols]
                        bond = mob[bond_index]
                        bond.put_start_and_end_on(atom.get_center(), neighbor.get_center())
                        bond.set_opacity(bond_opacity)
                        bond_index += 1
        
        # Updater for curve tracing - removed, curve is stationary
        
        # Updater for moving dot
        def update_moving_dot(mob):
            current_strain = strain.get_value()
            current_stress = stress_strain_curve(current_strain)
            mob.move_to(axes.coords_to_point(current_strain, current_stress))
        
        # Add updaters
        lattice_atoms.add_updater(update_lattice)
        bonds.add_updater(update_bonds)
        moving_dot.add_updater(update_moving_dot)
        
        # Animate the stretching
        self.play(
            strain.animate.set_value(max_strain),
            run_time=total_animation_time,
            rate_func=linear
        )
        
        # Remove updaters
        lattice_atoms.remove_updater(update_lattice)
        bonds.remove_updater(update_bonds)
        moving_dot.remove_updater(update_moving_dot)
        
        # Mark fracture point on graph
        
        # Split the lattice at fracture
        left_atoms = VGroup()
        right_atoms = VGroup()
        
        for i in range(rows):
            for j in range(cols):
                atom_index = i * cols + j
                if j < cols // 2:
                    left_atoms.add(lattice_atoms[atom_index])
                else:
                    right_atoms.add(lattice_atoms[atom_index])
        
        # Fade out bonds and separate halves
        self.play(
            FadeOut(bonds),
            left_atoms.animate.shift(0.5*LEFT),
            right_atoms.animate.shift(0.5*RIGHT),
            run_time=1.5,
            rate_func=rush_from
        )
        
        self.wait(2)
