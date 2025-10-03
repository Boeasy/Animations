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
        # divider = Line(UP, DOWN, color=WHITE, stroke_width=2)
        # divider.shift(0.3*RIGHT)
        
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
        ).shift(3.0*RIGHT + 0.3*DOWN)
        
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
            #Create(divider),
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
        
        # Animate the stretching up to fracture point
        time_to_fracture = total_animation_time * (fracture_strain / max_strain)
        self.play(
            strain.animate.set_value(fracture_strain),
            run_time=time_to_fracture,
            rate_func=linear
        )
        
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
        
        # Fade out bonds and separate halves (this happens AT the fracture point)
        self.play(
            FadeOut(bonds),
            left_atoms.animate.shift(0.5*LEFT),
            right_atoms.animate.shift(0.5*RIGHT),
            run_time=1.5,
            rate_func=rush_from
        )
        
        # Continue moving dot down the stress drop
        time_after_fracture = total_animation_time * ((max_strain - fracture_strain) / max_strain)
        self.play(
            strain.animate.set_value(max_strain),
            run_time=time_after_fracture,
            rate_func=linear
        )
        
        # Remove updaters
        lattice_atoms.remove_updater(update_lattice)
        bonds.remove_updater(update_bonds)
        moving_dot.remove_updater(update_moving_dot)
        
        self.wait(2)


class GoldBarFracture(Scene):
    def construct(self):
        # Create title
        title = Text("Stress-Strain Behavior: Gold Bar Fracture", font_size=32)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Create split screen divider
        # divider = Line(ORIGIN, 4*DOWN, color=WHITE, stroke_width=2)
        # divider.shift(0.3*RIGHT)
        
        # Left side: Gold Bar
        # Right side: Graph
        
        # Create macroscopic gold bar (trapezoidal prism like a gold bar)
        # Define vertices for a trapezoidal prism - oriented horizontally (long way)
        # Left trapezoid end (larger)
        left_front_bottom = [-2.0, -0.5, 0]
        left_front_top = [-2.0, 0.5, 0]
        left_back_bottom = [-2.0, -0.4, -0.4]
        left_back_top = [-2.0, 0.4, -0.4]
        
        # Right trapezoid end (smaller)
        right_front_bottom = [2.0, -0.4, 0]
        right_front_top = [2.0, 0.4, 0]
        right_back_bottom = [2.0, -0.3, -0.3]
        right_back_top = [2.0, 0.3, -0.3]
        
        vertices = [
            left_front_bottom, left_front_top, left_back_top, left_back_bottom,  # Left face
            right_front_bottom, right_front_top, right_back_top, right_back_bottom  # Right face
        ]
        
        # Define faces (each face is a list of vertex indices)
        faces = [
            [0, 1, 2, 3],  # Left trapezoid
            [4, 5, 6, 7],  # Right trapezoid
            [0, 1, 5, 4],  # Front face (bottom)
            [1, 2, 6, 5],  # Top face
            [2, 3, 7, 6],  # Back face
            [3, 0, 4, 7],  # Bottom face
        ]
        
        gold_bar = Polyhedron(
            vertex_coords=vertices,
            faces_list=faces,
            graph_config={"edge_config": {"color": GOLD_D, "stroke_width": 3}}
        )
        gold_bar.set_fill(GOLD, opacity=0.7)
        gold_bar.set_stroke(GOLD_D, width=2)
        
        # Rotate to get an angled view looking down at the gold bar
        gold_bar.rotate(angle=25*DEGREES, axis=RIGHT)  # Tilt forward
        gold_bar.rotate(angle=20*DEGREES, axis=UP)     # Rotate slightly for perspective
        
        # Shift gold bar to left side
        gold_bar.shift(3*LEFT + 0.3*DOWN)
        
        # Set up axes for stress-strain curve
        axes = Axes(
            x_range=[0, 0.25, 0.05],
            y_range=[0, 1.2, 0.2],
            x_length=4,
            y_length=3.5,
            axis_config={"include_tip": True, "numbers_to_exclude": []},
            tips=True,
        ).shift(3*RIGHT + 0.3*DOWN)
        
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
            #Create(divider),
            Create(gold_bar),
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
        
        # Store original vertices for deformation
        original_vertices = [np.array(v) for v in vertices]
        
        # Updater for gold bar stretching
        def update_gold_bar(mob):
            current_strain = strain.get_value()
            
            # Apply horizontal strain to each vertex (stretch along x-axis, the length)
            new_vertices = []
            for i, orig_vertex in enumerate(original_vertices):
                x_orig, y_orig, z_orig = orig_vertex
                
                # Apply horizontal strain (stretch in x-direction along the length)
                x_new = x_orig * (1 + current_strain * 2)  # Exaggerate for visibility
                
                # Add some vertical/depth contraction (Poisson effect)
                poisson_ratio = 0.44  # Gold's Poisson ratio
                y_new = y_orig * (1 - poisson_ratio * current_strain * 0.5)
                z_new = z_orig * (1 - poisson_ratio * current_strain * 0.5)
                
                new_vertices.append([x_new, y_new, z_new])
            
            # Recreate the polyhedron with new vertices
            new_bar = Polyhedron(
                vertex_coords=new_vertices,
                faces_list=faces,
                graph_config={"edge_config": {"color": GOLD_D, "stroke_width": 3}}
            )
            new_bar.set_fill(GOLD, opacity=0.7)
            new_bar.set_stroke(GOLD_D, width=2)
            new_bar.rotate(angle=25*DEGREES, axis=RIGHT)
            new_bar.rotate(angle=20*DEGREES, axis=UP)
            new_bar.shift(3*LEFT + 0.3*DOWN)
            
            # Fade effect as approaching fracture
            if current_strain > fracture_strain:
                opacity = max(0.3, 0.7 - 2*(current_strain - fracture_strain))
                new_bar.set_fill(opacity=opacity)
            
            mob.become(new_bar)
        
        # Updater for moving dot
        def update_moving_dot(mob):
            current_strain = strain.get_value()
            current_stress = stress_strain_curve(current_strain)
            mob.move_to(axes.coords_to_point(current_strain, current_stress))
        
        # Add updaters
        gold_bar.add_updater(update_gold_bar)
        moving_dot.add_updater(update_moving_dot)
        
        # Animate the stretching up to fracture point
        time_to_fracture = total_animation_time * (fracture_strain / max_strain)
        self.play(
            strain.animate.set_value(fracture_strain),
            run_time=time_to_fracture,
            rate_func=linear
        )
        
        # Create two halves of the fractured bar (split vertically down the middle)
        # Use fracture_strain values for the split
        poisson_y = (1 - 0.44 * fracture_strain * 0.5)
        poisson_z = (1 - 0.44 * fracture_strain * 0.5)
        stretch = (1 + fracture_strain * 2)
        
        left_vertices = [
            [-2.0 * stretch, -0.5 * poisson_y, 0],
            [-2.0 * stretch, 0.5 * poisson_y, 0],
            [-2.0 * stretch, 0.4 * poisson_y, -0.4 * poisson_z],
            [-2.0 * stretch, -0.4 * poisson_y, -0.4 * poisson_z],
            [0, -0.4 * poisson_y, 0],
            [0, 0.4 * poisson_y, 0],
            [0, 0.3 * poisson_y, -0.3 * poisson_z],
            [0, -0.3 * poisson_y, -0.3 * poisson_z],
        ]
        
        # Right half
        right_vertices = [
            [0, -0.4 * poisson_y, 0],
            [0, 0.4 * poisson_y, 0],
            [0, 0.3 * poisson_y, -0.3 * poisson_z],
            [0, -0.3 * poisson_y, -0.3 * poisson_z],
            [2.0 * stretch, -0.4 * poisson_y, 0],
            [2.0 * stretch, 0.4 * poisson_y, 0],
            [2.0 * stretch, 0.3 * poisson_y, -0.3 * poisson_z],
            [2.0 * stretch, -0.3 * poisson_y, -0.3 * poisson_z],
        ]
        
        left_half = Polyhedron(
            vertex_coords=left_vertices,
            faces_list=faces,
            graph_config={"edge_config": {"color": GOLD_D, "stroke_width": 3}}
        )
        left_half.set_fill(GOLD, opacity=0.7)
        left_half.set_stroke(GOLD_D, width=2)
        left_half.rotate(angle=25*DEGREES, axis=RIGHT)
        left_half.rotate(angle=20*DEGREES, axis=UP)
        left_half.shift(3*LEFT + 0.3*DOWN)
        
        right_half = Polyhedron(
            vertex_coords=right_vertices,
            faces_list=faces,
            graph_config={"edge_config": {"color": GOLD_D, "stroke_width": 3}}
        )
        right_half.set_fill(GOLD, opacity=0.7)
        right_half.set_stroke(GOLD_D, width=2)
        right_half.rotate(angle=25*DEGREES, axis=RIGHT)
        right_half.rotate(angle=20*DEGREES, axis=UP)
        right_half.shift(3*LEFT + 0.3*DOWN)
        
        # Replace the stretched bar with the two halves and separate them (this happens AT fracture point)
        self.play(
            FadeOut(gold_bar),
            FadeIn(left_half),
            FadeIn(right_half),
            run_time=0.3
        )
        
        self.play(
            left_half.animate.shift(0.8*LEFT),
            right_half.animate.shift(0.8*RIGHT),
            run_time=1.5,
            rate_func=rush_from
        )
        
        # Continue moving dot down the stress drop
        time_after_fracture = total_animation_time * ((max_strain - fracture_strain) / max_strain)
        self.play(
            strain.animate.set_value(max_strain),
            run_time=time_after_fracture,
            rate_func=linear
        )
        
        # Remove updaters
        gold_bar.remove_updater(update_gold_bar)
        moving_dot.remove_updater(update_moving_dot)
        
        self.wait(2)
