""" Visualizing the stress strain curve of gold (or a general metal) up to the fracture while stretching, 
visualizing the graph and the atomistic behavior

In good shape, can't get the force vectors to work with this, keeps having weird errors (not updating unless reboot)
 """

from manim import *
import numpy as np

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
                    
                    # Check if this is a center column bond (between cols//2-1 and cols//2)
                    is_center_bond = (j == cols // 2 - 1)
                    
                    # Original position
                    x_orig = (j - cols/2 + 0.5) * initial_spacing
                    y_orig = (i - rows/2 + 0.5) * initial_spacing - 0.5
                    
                    # Apply horizontal strain from center
                    distance_from_center = x_orig
                    if is_center_bond:
                        #Exaggerate the center bond stretch 
                        x_new = x_orig * (1 + 10*current_strain)
                    else:
                        x_new = x_orig * (1 + current_strain)
                    
                    # Add some vertical contraction (Poisson effect)
                    poisson_ratio = 0.44  # Gold's Poisson ratio
                    y_new = y_orig * (1 - poisson_ratio * current_strain * 0.5)
                    
                    atom.move_to([x_new - 3, y_new, 0])
        
        # Updater for bonds
        def update_bonds(mob):
            """ This function will update the bonds of each of the atoms in the lattice.
            
            We need to loop over all atoms and update the bonds each time, because the atoms are moving from the lattice updater above.

            When the fracture point is reached, I want the bonds to sequentially break, top to bottom. 

            """
            current_strain = strain.get_value()
            bond_index = 0
            
            for i in range(rows):
                for j in range(cols):
                    atom_index = i * cols + j
                    atom = lattice_atoms[atom_index]
                    
                    # Horizontal bonds
                    if j < cols - 1:
                        neighbor = lattice_atoms[atom_index + 1]
                        bond = mob[bond_index]
                        
                        # Update bond position
                        bond.put_start_and_end_on(atom.get_center(), neighbor.get_center())
                        
                        # Check if this is a center column bond (between cols//2-1 and cols//2)
                        is_center_bond = (j == cols // 2 - 1)
                        
                        if is_center_bond:
                            # Change stroke width based on strain (thinner as it stretches)
                            base_width = 2
                            new_width = max(0.5, base_width * (1 - current_strain * 3))
                            bond.set_stroke(width=new_width)
                            
                            # Break bonds sequentially top to bottom when fracture is reached
                            if current_strain >= fracture_strain:
                                # Calculate which row should break based on how far past fracture we are
                                progress_past_fracture = (current_strain - fracture_strain) / (max_strain - fracture_strain)
                                rows_to_break = int(progress_past_fracture * rows)
                                
                                # Break from top (highest index) to bottom
                                row_from_top = rows - 1 - i
                                if row_from_top < rows_to_break:
                                    bond.set_opacity(0)
                        
                        bond_index += 1
                    
                    # Vertical bonds
                    if i < rows - 1:
                        neighbor = lattice_atoms[atom_index + cols]
                        bond = mob[bond_index]
                        
                        # Update bond position
                        bond.put_start_and_end_on(atom.get_center(), neighbor.get_center())
                        
                        bond_index += 1
                
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
        time_after_fracture = total_animation_time * ((max_strain - fracture_strain) / max_strain)
        self.play(
            strain.animate.set_value(max_strain),
            run_time=(time_to_fracture+time_after_fracture),
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
        
        self.play(
            FadeOut(bonds),
            left_atoms.animate.shift(0.3*LEFT),
            right_atoms.animate.shift(0.3*RIGHT),
            run_time=0.5,
            rate_func=rush_from
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
        
        a = 70
        b = 50

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
        gold_bar.rotate(angle=a*DEGREES, axis=RIGHT)  # Tilt forward
        gold_bar.rotate(angle=b*DEGREES, axis=UP)     # Rotate slightly for perspective
        
        # Shift gold bar to left side
        gold_bar.shift(3.5*LEFT + 0.3*DOWN)
        
        # Set up axes for stress-strain curve
        axes = Axes(
            x_range=[0, 0.25, 0.05],
            y_range=[0, 1.2, 0.2],
            x_length=4,
            y_length=3.5,
            axis_config={"include_tip": True, "numbers_to_exclude": []},
            tips=True,
        ).shift(3.5*RIGHT + 0.3*DOWN)
        
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

        def update_gold_bar(mob):
            current_strain = strain.get_value()
            
            # If we haven't split yet and strain is below fracture, show unified bar
            if current_strain < fracture_strain:
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
                new_bar.rotate(angle=a*DEGREES, axis=RIGHT)
                new_bar.rotate(angle=b*DEGREES, axis=UP)
                new_bar.shift(3.5*LEFT + 0.3*DOWN)
                
                mob.become(new_bar)
            
            else:
                # We're at or past fracture - create/update split halves
                # Calculate stretched vertices at fracture strain for the split
                stretched_vertices = []
                for orig_vertex in original_vertices:
                    x_orig, y_orig, z_orig = orig_vertex
                    x_new = x_orig * (1 + fracture_strain * 2)
                    poisson_ratio = 0.44
                    y_new = y_orig * (1 - poisson_ratio * fracture_strain * 0.5)
                    z_new = z_orig * (1 - poisson_ratio * fracture_strain * 0.5)
                    stretched_vertices.append(np.array([x_new, y_new, z_new]))
                
                # Get the stretched x-coordinates for interpolation
                left_x = stretched_vertices[0][0]
                right_x = stretched_vertices[4][0]
                
                # Interpolate center vertices at x=0
                center_front_bottom = np.array([0, 
                    np.interp(0, [left_x, right_x], [stretched_vertices[0][1], stretched_vertices[4][1]]),
                    np.interp(0, [left_x, right_x], [stretched_vertices[0][2], stretched_vertices[4][2]])])
                
                center_front_top = np.array([0,
                    np.interp(0, [left_x, right_x], [stretched_vertices[1][1], stretched_vertices[5][1]]),
                    np.interp(0, [left_x, right_x], [stretched_vertices[1][2], stretched_vertices[5][2]])])
                
                center_back_top = np.array([0,
                    np.interp(0, [left_x, right_x], [stretched_vertices[2][1], stretched_vertices[6][1]]),
                    np.interp(0, [left_x, right_x], [stretched_vertices[2][2], stretched_vertices[6][2]])])
                
                center_back_bottom = np.array([0,
                    np.interp(0, [left_x, right_x], [stretched_vertices[3][1], stretched_vertices[7][1]]),
                    np.interp(0, [left_x, right_x], [stretched_vertices[3][2], stretched_vertices[7][2]])])
                
                # Calculate separation based on how far past fracture we are
                separation_progress = (current_strain - fracture_strain) / (max_strain - fracture_strain)
                separation_distance = separation_progress * 0.9  # Total separation of 1.5 units
                
                # Build left half vertices with separation
                left_half_vertices = [
                    stretched_vertices[0] + np.array([-separation_distance, 0, 0]),
                    stretched_vertices[1] + np.array([-separation_distance, 0, 0]),
                    stretched_vertices[2] + np.array([-separation_distance, 0, 0]),
                    stretched_vertices[3] + np.array([-separation_distance, 0, 0]),
                    center_front_bottom + np.array([-separation_distance, 0, 0]),
                    center_front_top + np.array([-separation_distance, 0, 0]),
                    center_back_top + np.array([-separation_distance, 0, 0]),
                    center_back_bottom + np.array([-separation_distance, 0, 0]),
                ]
                
                # Build right half vertices with separation
                right_half_vertices = [
                    center_front_bottom + np.array([separation_distance, 0, 0]),
                    center_front_top + np.array([separation_distance, 0, 0]),
                    center_back_top + np.array([separation_distance, 0, 0]),
                    center_back_bottom + np.array([separation_distance, 0, 0]),
                    stretched_vertices[4] + np.array([separation_distance, 0, 0]),
                    stretched_vertices[5] + np.array([separation_distance, 0, 0]),
                    stretched_vertices[6] + np.array([separation_distance, 0, 0]),
                    stretched_vertices[7] + np.array([separation_distance, 0, 0]),
                ]
                
                # Create both halves
                left_half = Polyhedron(
                    vertex_coords=left_half_vertices,
                    faces_list=faces,
                    graph_config={"edge_config": {"color": GOLD_D, "stroke_width": 3}}
                )
                left_half.set_fill(GOLD, opacity=0.7)
                left_half.set_stroke(GOLD_D, width=2)
                left_half.rotate(angle=a*DEGREES, axis=RIGHT)
                left_half.rotate(angle=b*DEGREES, axis=UP)
                left_half.shift(3.5*LEFT + 0.3*DOWN)
                
                right_half = Polyhedron(
                    vertex_coords=right_half_vertices,
                    faces_list=faces,
                    graph_config={"edge_config": {"color": GOLD_D, "stroke_width": 3}}
                )
                right_half.set_fill(GOLD, opacity=0.7)
                right_half.set_stroke(GOLD_D, width=2)
                right_half.rotate(angle=a*DEGREES, axis=RIGHT)
                right_half.rotate(angle=b*DEGREES, axis=UP)
                right_half.shift(3.5*LEFT + 0.3*DOWN)
                
                # Combine both halves into a single VGroup
                split_bars = VGroup(left_half, right_half)
                mob.become(split_bars)
        
        # Updater for moving dot
        def update_moving_dot(mob):
            current_strain = strain.get_value()
            current_stress = stress_strain_curve(current_strain)
            mob.move_to(axes.coords_to_point(current_strain, current_stress))
        
        # Add updaters
        gold_bar.add_updater(update_gold_bar)
        moving_dot.add_updater(update_moving_dot)
        
        # Animate the entire stretching and fracture in one continuous motion
        self.play(
            strain.animate.set_value(max_strain),
            run_time=total_animation_time,
            rate_func=linear
        )
        
        # Remove updaters
        gold_bar.remove_updater(update_gold_bar)
        moving_dot.remove_updater(update_moving_dot)
        
        self.wait(2)
