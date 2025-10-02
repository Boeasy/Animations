"""This scene will frame how molecular dynamics help us find macroscopic properties of materials"""
#Pretty close, fix up the gold prism, move the zoom to something that is actually a part of the prism lol, and then make sure the atoms don't move around so much in the lattice. maybe leave the lattice stationary and let the atoms move around around those points.


from manim import *
import numpy as np


class MDMacro(Scene):
    def construct(self):
        # Title
        title = Text("From Macro to Micro: Understanding Temperature", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Create macroscopic gold bar (rectangular prism in isometric view)
        bar_width = 4.0
        bar_height = 2.5
        bar_depth = 1.5
        
        def apply_iso(x, y, z):
            """Apply isometric projection"""
            iso_x = x - z * 0.6
            iso_y = y + (x + z) * 0.3
            return np.array([iso_x, iso_y, 0])
        
        # Create vertices of the rectangular prism
        vertices_3d = [
            np.array([0, 0, 0]),
            np.array([bar_width, 0, 0]),
            np.array([bar_width, bar_height, 0]),
            np.array([0, bar_height, 0]),
            np.array([0, 0, bar_depth]),
            np.array([bar_width, 0, bar_depth]),
            np.array([bar_width, bar_height, bar_depth]),
            np.array([0, bar_height, bar_depth]),
        ]
        
        # Apply isometric projection and center
        center_3d = np.array([bar_width/2, bar_height/2, bar_depth/2])
        vertices_2d = [apply_iso(*(v - center_3d)) for v in vertices_3d]
        
        # Create faces
        front_face = Polygon(*[vertices_2d[i] for i in [0, 1, 2, 3]], 
                            color=GOLD, fill_opacity=0.8, stroke_color=YELLOW, stroke_width=2)
        right_face = Polygon(*[vertices_2d[i] for i in [1, 5, 6, 2]], 
                            color=GOLD, fill_opacity=0.6, stroke_color=YELLOW, stroke_width=2)
        top_face = Polygon(*[vertices_2d[i] for i in [3, 2, 6, 7]], 
                          color=GOLD, fill_opacity=0.9, stroke_color=YELLOW, stroke_width=2)
        
        gold_bar = VGroup(front_face, right_face, top_face)
        gold_bar.shift(0.5*DOWN)
        
        # Add label
        bar_label = Text("Gold Bar", font_size=28, color=YELLOW, weight=BOLD)
        bar_label.next_to(gold_bar, DOWN, buff=0.5)
        
        # Show the gold bar
        self.play(
            Create(gold_bar),
            Write(bar_label),
            run_time=2
        )
        self.wait(0.5)
        
        # Add temperature label
        temp_value = 300  # Kelvin
        temp_label = VGroup(
            Text("Temperature:", font_size=24),
            MathTex("T = 300", r"\,\text{K}", font_size=28, color=RED)
        ).arrange(RIGHT, buff=0.3)
        temp_label.next_to(gold_bar, LEFT, buff=1)
        
        self.play(Write(temp_label), run_time=1)
        self.wait(0.5)
        
        # Add other observables
        density_label = VGroup(
            Text("Density:", font_size=20),
            MathTex(r"\rho = 19.3\,\text{g/cm}^3", font_size=22, color=BLUE)
        ).arrange(RIGHT, buff=0.2)
        density_label.next_to(temp_label, DOWN, aligned_edge=LEFT, buff=0.3)
        
        mass_label = VGroup(
            Text("Mass:", font_size=20),
            MathTex("m = 150", r"\,\text{g}", font_size=22, color=GREEN)
        ).arrange(RIGHT, buff=0.2)
        mass_label.next_to(density_label, DOWN, aligned_edge=LEFT, buff=0.3)
        
        self.play(
            Write(density_label),
            Write(mass_label),
            run_time=1.5
        )
        self.wait(1)
        
        # Highlight a region to zoom into
        zoom_region = Rectangle(
            width=1.0, 
            height=0.8,
            color=WHITE,
            stroke_width=4
        )
        zoom_region.move_to(apply_iso(0.3, 1.0, 0.3))
        
        zoom_label = Text("Zoom in...", font_size=20, color=WHITE)
        zoom_label.next_to(zoom_region, RIGHT, buff=0.3)
        
        self.play(
            Create(zoom_region),
            Write(zoom_label),
            run_time=1
        )
        self.wait(0.5)
        
        # Fade out everything except zoom region
        self.play(
            FadeOut(bar_label),
            FadeOut(temp_label),
            FadeOut(density_label),
            FadeOut(mass_label),
            FadeOut(zoom_label),
            gold_bar.animate.set_opacity(0.3),
            run_time=1
        )
        self.wait(0.3)
        
        # Create atomic lattice structure
        lattice_size = 5
        atom_spacing = 0.6
        atom_radius = 0.12
        
        # Create FCC-like arrangement (simplified to 2D view)
        atoms = VGroup()
        atom_positions = []
        
        for i in range(lattice_size):
            for j in range(lattice_size):
                x = (i - lattice_size/2 + 0.5) * atom_spacing
                y = (j - lattice_size/2 + 0.5) * atom_spacing
                
                atom = Circle(
                    radius=atom_radius,
                    color=GOLD,
                    fill_opacity=1,
                    stroke_width=2,
                    stroke_color=YELLOW
                )
                atom.move_to([x, y, 0])
                atoms.add(atom)
                atom_positions.append(np.array([x, y, 0]))
        
        # Add bonds between atoms
        bonds = VGroup()
        for i in range(lattice_size):
            for j in range(lattice_size):
                idx = i * lattice_size + j
                
                # Horizontal bonds
                if i < lattice_size - 1:
                    neighbor_idx = (i + 1) * lattice_size + j
                    bond = Line(
                        atom_positions[idx],
                        atom_positions[neighbor_idx],
                        color=GRAY,
                        stroke_width=2,
                        stroke_opacity=0.5
                    )
                    bonds.add(bond)
                
                # Vertical bonds
                if j < lattice_size - 1:
                    neighbor_idx = i * lattice_size + (j + 1)
                    bond = Line(
                        atom_positions[idx],
                        atom_positions[neighbor_idx],
                        color=GRAY,
                        stroke_width=2,
                        stroke_opacity=0.5
                    )
                    bonds.add(bond)
        
        # Zoom transition
        lattice_group = VGroup(bonds, atoms)
        lattice_group.scale(0.1).move_to(zoom_region.get_center())
        
        self.play(
            FadeOut(gold_bar),
            FadeOut(zoom_region),
            lattice_group.animate.scale(10).move_to(ORIGIN).shift(0.5*UP),
            run_time=2
        )
        
        # Update title
        self.play(
            Transform(title, Text("Atomic Scale: Kinetic Energy & Temperature", font_size=32).to_edge(UP)),
            run_time=1
        )
        self.wait(0.5)
        
        # Add thermal vibration to atoms
        np.random.seed(42)
        vibration_amplitude = 0.08
        
        # Store original positions
        original_positions = [atom.get_center().copy() for atom in atoms]
        
        # Create velocity vectors for visualization (just a few to avoid clutter)
        show_indices = [6, 8, 12, 16, 18]  # Show vectors for select atoms
        velocity_arrows = VGroup()
        
        for idx in show_indices:
            arrow = Arrow(
                ORIGIN,
                RIGHT * 0.5,
                color=RED,
                buff=0,
                stroke_width=3,
                max_tip_length_to_length_ratio=0.2
            )
            velocity_arrows.add(arrow)
        
        # Function to update atom positions with thermal motion
        time_tracker = ValueTracker(0)
        
        def update_atoms(mob):
            t = time_tracker.get_value()
            for i, atom in enumerate(mob):
                # Each atom vibrates with its own phase
                phase_x = i * 0.7
                phase_y = i * 1.3
                
                offset_x = vibration_amplitude * np.sin(3 * t + phase_x)
                offset_y = vibration_amplitude * np.sin(3 * t + phase_y)
                
                new_pos = original_positions[i] + np.array([offset_x, offset_y, 0])
                atom.move_to(new_pos)
        
        def update_bonds(mob):
            bond_idx = 0
            for i in range(lattice_size):
                for j in range(lattice_size):
                    idx = i * lattice_size + j
                    
                    # Horizontal bonds
                    if i < lattice_size - 1:
                        neighbor_idx = (i + 1) * lattice_size + j
                        mob[bond_idx].put_start_and_end_on(
                            atoms[idx].get_center(),
                            atoms[neighbor_idx].get_center()
                        )
                        bond_idx += 1
                    
                    # Vertical bonds
                    if j < lattice_size - 1:
                        neighbor_idx = i * lattice_size + (j + 1)
                        mob[bond_idx].put_start_and_end_on(
                            atoms[idx].get_center(),
                            atoms[neighbor_idx].get_center()
                        )
                        bond_idx += 1
        
        def update_velocity_arrows(mob):
            t = time_tracker.get_value()
            for arrow_idx, atom_idx in enumerate(show_indices):
                atom = atoms[atom_idx]
                
                # Calculate instantaneous velocity direction
                phase_x = atom_idx * 0.7
                phase_y = atom_idx * 1.3
                
                vel_x = vibration_amplitude * 3 * np.cos(3 * t + phase_x)
                vel_y = vibration_amplitude * 3 * np.cos(3 * t + phase_y)
                
                velocity = np.array([vel_x, vel_y, 0])
                speed = np.linalg.norm(velocity)
                
                if speed > 0.01:
                    direction = velocity / speed
                    arrow_length = min(speed * 2, 0.6)
                    
                    mob[arrow_idx].put_start_and_end_on(
                        atom.get_center(),
                        atom.get_center() + direction * arrow_length
                    )
                    mob[arrow_idx].set_opacity(1)
                else:
                    mob[arrow_idx].set_opacity(0)
        
        # Add updaters
        atoms.add_updater(update_atoms)
        bonds.add_updater(update_bonds)
        
        # Show atoms vibrating
        self.add(velocity_arrows)
        velocity_arrows.add_updater(update_velocity_arrows)
        
        velocity_label = Text("Velocity vectors", font_size=20, color=RED)
        velocity_label.to_edge(RIGHT).shift(UP)
        
        self.play(
            Write(velocity_label),
            run_time=1
        )
        
        # Animate thermal motion
        self.play(
            time_tracker.animate.set_value(4),
            run_time=4,
            rate_func=linear
        )
        
        # Show kinetic energy formula
        ke_formula = MathTex(
            r"\langle KE \rangle = \frac{1}{2}m\langle v^2 \rangle",
            font_size=32,
            color=YELLOW
        )
        ke_formula.to_edge(LEFT).shift(UP)
        
        self.play(Write(ke_formula), run_time=1)
        
        # Continue motion
        self.play(
            time_tracker.animate.set_value(8),
            run_time=4,
            rate_func=linear
        )
        
        # Show connection to temperature
        connection_arrow = Arrow(
            ke_formula.get_bottom(),
            ke_formula.get_bottom() + DOWN * 0.8,
            color=WHITE,
            stroke_width=4
        )
        
        temp_connection = MathTex(
            r"\langle KE \rangle = \frac{3}{2}k_B T",
            font_size=32,
            color=RED
        )
        temp_connection.next_to(connection_arrow, DOWN)
        
        self.play(
            GrowArrow(connection_arrow),
            Write(temp_connection),
            run_time=1.5
        )
        
        # Continue motion
        self.play(
            time_tracker.animate.set_value(10),
            run_time=2,
            rate_func=linear
        )
        
        # Final explanation
        explanation = VGroup(
            Text("Macroscopic Temperature", font_size=24, color=RED, weight=BOLD),
            MathTex(r"\Updownarrow", font_size=32),
            Text("Average Kinetic Energy", font_size=24, color=YELLOW, weight=BOLD),
            Text("of microscopic atoms", font_size=20, color=YELLOW)
        ).arrange(DOWN, buff=0.2)
        explanation.to_edge(RIGHT).shift(0.5*DOWN)
        
        self.play(Write(explanation), run_time=2)
        
        # Continue vibration
        self.play(
            time_tracker.animate.set_value(14),
            run_time=4,
            rate_func=linear
        )
        
        self.wait(2)