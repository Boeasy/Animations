from manim import *
import numpy as np


class MLPotential(Scene):
    """Demonstrate the difference between EAM and ML potentials."""
    
    def construct(self):
        # Title
        title = Text("ML Potential vs. EAM", font_size=40)
        title.to_edge(UP)
        self.play(FadeIn(title))
        self.wait(0.5)
        
        # Parameters
        lattice_size = 5  # 5x5 lattice by default
        atom_radius = 0.18
        initial_spacing = 1.5  # Start further apart
        final_spacing = 1.0    # Final lattice spacing
        cloud_radius = 0.45
        gold_color = GOLD
        cloud_color = YELLOW
        pair_color = GREEN
        repulsion_color = RED
        
        # ============================================================
        # PART 1: PAIR Potential - creates lattice
        # ============================================================
        """
        We're going to create a 5x5 lattice here of gold atoms. They'll start somewhat further apart, and move in to a lattice structure.
        """
        
        # Label for Part 1
        part1_label = Text("Pair Potential", font_size=28, color=pair_color)
        part1_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(part1_label))
        self.wait(0.3)
        
        # Create lattice positions (starting far apart)
        atoms_part1 = VGroup()
        initial_positions = []
        final_positions = []
        
        center_offset = np.array([0, -0.5, 0])  # Shift lattice down slightly
        
        for i in range(lattice_size):
            for j in range(lattice_size):
                # Initial position (spread out)
                x_init = (j - lattice_size/2 + 0.5) * initial_spacing
                y_init = (i - lattice_size/2 + 0.5) * initial_spacing
                initial_pos = np.array([x_init, y_init, 0]) + center_offset
                initial_positions.append(initial_pos)
                
                # Final position (lattice)
                x_final = (j - lattice_size/2 + 0.5) * final_spacing
                y_final = (i - lattice_size/2 + 0.5) * final_spacing
                final_pos = np.array([x_final, y_final, 0]) + center_offset
                final_positions.append(final_pos)
                
                # Create atom
                atom = Circle(radius=atom_radius, stroke_width=3)
                atom.set_fill(gold_color, opacity=1.0)
                atom.set_stroke(GOLD_E, width=3)
                atom.move_to(initial_pos)
                atoms_part1.add(atom)
        
        # Show atoms appearing
        self.play(
            LaggedStart(*[FadeIn(atom, scale=0.5) for atom in atoms_part1], lag_ratio=0.03),
            run_time=2.0
        )
        self.wait(0.5)
        
        # Draw lattice bonds between atoms at initial positions
        def create_lattice_bonds(positions, spacing_val):
            bonds = VGroup()
            for i in range(lattice_size):
                for j in range(lattice_size):
                    idx = i * lattice_size + j
                    # Horizontal bonds
                    if j < lattice_size - 1:
                        start = positions[idx]
                        end = positions[idx + 1]
                        bond = Line(start, end, color=WHITE, stroke_width=2)
                        bonds.add(bond)
                    # Vertical bonds
                    if i < lattice_size - 1:
                        start = positions[idx]
                        end = positions[idx + lattice_size]
                        bond = Line(start, end, color=WHITE, stroke_width=2)
                        bonds.add(bond)
            return bonds
        
        # Move atoms together into lattice
        animations = []
        for i, atom in enumerate(atoms_part1):
            animations.append(atom.animate.move_to(final_positions[i]))
        
        self.play(*animations, run_time=2.5, rate_func=smooth)
        self.wait(0.5)
        
        # Draw lattice bonds
        lattice_bonds_part1 = create_lattice_bonds(final_positions, final_spacing)
        self.play(Create(lattice_bonds_part1), run_time=1.5)
        self.wait(1.0)
        
        # Fade out everything from Part 1
        self.play(
            FadeOut(atoms_part1),
            FadeOut(lattice_bonds_part1),
            FadeOut(part1_label),
            run_time=1.0
        )
        self.wait(0.5)
        
        # ============================================================
        # PART 2: EAM Potential with symmetric electron clouds
        # ============================================================
        """
        EAM Potential adds the circle electron clouds to the lattice, and then the repulsive force of the electron clouds moves the lattice apart slightly.
        
        electron clouds fade out at end.
        """
        
        # Label for Part 2
        part2_label = Text("EAM Potential", font_size=28, color=cloud_color)
        part2_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(part2_label))
        self.wait(0.3)
        
        # Create new atoms at final lattice positions
        atoms_part2 = VGroup()
        for pos in final_positions:
            atom = Circle(radius=atom_radius, stroke_width=3)
            atom.set_fill(gold_color, opacity=1.0)
            atom.set_stroke(GOLD_E, width=3)
            atom.move_to(pos)
            atoms_part2.add(atom)
        
        self.play(FadeIn(atoms_part2, scale=0.8), run_time=1.0)
        self.wait(0.3)
        
        # Create symmetric electron clouds
        electron_clouds = VGroup()
        for pos in final_positions:
            cloud = Circle(radius=cloud_radius, stroke_width=0)
            cloud.set_fill(cloud_color, opacity=0.2)
            cloud.set_stroke(cloud_color, width=1, opacity=0.4)
            cloud.move_to(pos)
            electron_clouds.add(cloud)
        
        # Show electron clouds
        self.play(
            LaggedStart(*[FadeIn(cloud, scale=0.5) for cloud in electron_clouds], lag_ratio=0.03),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Move atoms apart slightly due to electron repulsion
        eam_spacing = final_spacing * 1.15  # Slightly larger spacing
        eam_positions = []
        
        for i in range(lattice_size):
            for j in range(lattice_size):
                x_eam = (j - lattice_size/2 + 0.5) * eam_spacing
                y_eam = (i - lattice_size/2 + 0.5) * eam_spacing
                eam_pos = np.array([x_eam, y_eam, 0]) + center_offset
                eam_positions.append(eam_pos)
        
        # Animate atoms and clouds moving apart
        animations = []
        for i, (atom, cloud) in enumerate(zip(atoms_part2, electron_clouds)):
            animations.append(atom.animate.move_to(eam_positions[i]))
            animations.append(cloud.animate.move_to(eam_positions[i]))
        
        self.play(*animations, run_time=2.0, rate_func=smooth)
        self.wait(0.5)
        
        # Draw lattice bonds at EAM positions
        lattice_bonds_part2 = create_lattice_bonds(eam_positions, eam_spacing)
        self.play(Create(lattice_bonds_part2), run_time=1.5)
        self.wait(1.0)
        
        # Fade out everything from Part 2
        self.play(
            FadeOut(atoms_part2),
            FadeOut(electron_clouds),
            FadeOut(lattice_bonds_part2),
            FadeOut(part2_label),
            run_time=1.0
        )
        self.wait(0.5)
        
        # ============================================================
        # PART 3: ML Potential with non-symmetric electron clouds
        # ============================================================
        """
        ML Potential changes the electron clouds to s-d hybridization orbital shapes, changing the lattice in an anisotropic way.
        """
        
        # Label for Part 3
        part3_label = Text("Machine Learning Potential", font_size=28, color=RED)
        part3_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(part3_label))
        self.wait(0.3)
        
        # Create new atoms at final lattice positions
        atoms_part3 = VGroup()
        for pos in final_positions:
            atom = Circle(radius=atom_radius, stroke_width=3)
            atom.set_fill(gold_color, opacity=1.0)
            atom.set_stroke(GOLD_E, width=3)
            atom.move_to(pos)
            atoms_part3.add(atom)
        
        self.play(FadeIn(atoms_part3, scale=0.8), run_time=1.0)
        self.wait(0.3)
        
        # Create non-symmetric electron clouds (ellipses, wider than tall)
        hybrid_clouds = VGroup()
        for pos in final_positions:
            # Ellipse oriented horizontally (s-d hybrid orbitals)
            cloud = Ellipse(
                width=cloud_radius * 2.5,   # Wider
                height=cloud_radius * 1.3,  # Narrower
                stroke_width=0
            )
            cloud.set_fill(cloud_color, opacity=0.2)
            cloud.set_stroke(cloud_color, width=1, opacity=0.4)
            cloud.move_to(pos)
            hybrid_clouds.add(cloud)
        
        # Show hybrid electron clouds
        self.play(
            LaggedStart(*[FadeIn(cloud, scale=0.5) for cloud in hybrid_clouds], lag_ratio=0.03),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Move atoms in anisotropic way (more spread horizontally than vertically)
        ml_spacing_x = final_spacing * 1.3  # Spread more left-right
        ml_spacing_y = final_spacing * 1.05  # Spread less up-down
        ml_positions = []
        
        for i in range(lattice_size):
            for j in range(lattice_size):
                x_ml = (j - lattice_size/2 + 0.5) * ml_spacing_x
                y_ml = (i - lattice_size/2 + 0.5) * ml_spacing_y
                ml_pos = np.array([x_ml, y_ml, 0]) + center_offset
                ml_positions.append(ml_pos)
        
        # Animate atoms and clouds moving anisotropically
        animations = []
        for i, (atom, cloud) in enumerate(zip(atoms_part3, hybrid_clouds)):
            animations.append(atom.animate.move_to(ml_positions[i]))
            animations.append(cloud.animate.move_to(ml_positions[i]))
        
        self.play(*animations, run_time=2.0, rate_func=smooth)
        self.wait(0.5)
        
        # Draw lattice bonds at ML positions
        def create_anisotropic_lattice_bonds(positions):
            bonds = VGroup()
            for i in range(lattice_size):
                for j in range(lattice_size):
                    idx = i * lattice_size + j
                    # Horizontal bonds
                    if j < lattice_size - 1:
                        start = positions[idx]
                        end = positions[idx + 1]
                        bond = Line(start, end, color=WHITE, stroke_width=2)
                        bonds.add(bond)
                    # Vertical bonds
                    if i < lattice_size - 1:
                        start = positions[idx]
                        end = positions[idx + lattice_size]
                        bond = Line(start, end, color=WHITE, stroke_width=2)
                        bonds.add(bond)
            return bonds
        
        lattice_bonds_part3 = create_anisotropic_lattice_bonds(ml_positions)
        self.play(Create(lattice_bonds_part3), run_time=1.5)
        self.wait(1.0)
        
        # Fade out current display
        self.play(
            FadeOut(atoms_part3),
            FadeOut(hybrid_clouds),
            FadeOut(lattice_bonds_part3),
            FadeOut(part3_label),
            run_time=1.0
        )
        self.wait(0.5)
        
        # ============================================================
        # PART 4: Comparison
        # ============================================================
        """
        The three lattices at the end of their animations are now compared side by side.
        """
        
        # Label for comparison
        comparison_label = Text("Comparison", font_size=32, color=WHITE)
        comparison_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(comparison_label))
        self.wait(0.5)
        
        # Scale factor for side-by-side display
        scale_factor = 0.6
        
        # Create three lattices side by side
        # Part 1 lattice (Pair Potential) - LEFT
        atoms_final_1 = VGroup()
        bonds_final_1 = VGroup()
        for i, pos in enumerate(final_positions):
            atom = Circle(radius=atom_radius * scale_factor, stroke_width=2)
            atom.set_fill(gold_color, opacity=1.0)
            atom.set_stroke(GOLD_E, width=2)
            # Scale and shift to left
            scaled_pos = (pos - center_offset) * scale_factor + np.array([-3.5, -0.5, 0])
            atom.move_to(scaled_pos)
            atoms_final_1.add(atom)
        
        # Bonds for Part 1
        for i in range(lattice_size):
            for j in range(lattice_size):
                idx = i * lattice_size + j
                atom_pos = atoms_final_1[idx].get_center()
                if j < lattice_size - 1:
                    neighbor_pos = atoms_final_1[idx + 1].get_center()
                    bond = Line(atom_pos, neighbor_pos, color=WHITE, stroke_width=1.5)
                    bonds_final_1.add(bond)
                if i < lattice_size - 1:
                    neighbor_pos = atoms_final_1[idx + lattice_size].get_center()
                    bond = Line(atom_pos, neighbor_pos, color=WHITE, stroke_width=1.5)
                    bonds_final_1.add(bond)
        
        label_1 = Text("Pair", font_size=20, color=pair_color)
        label_1.next_to(atoms_final_1, UP, buff=0.3)
        
        # Part 2 lattice (EAM) - CENTER
        atoms_final_2 = VGroup()
        bonds_final_2 = VGroup()
        for i, pos in enumerate(eam_positions):
            atom = Circle(radius=atom_radius * scale_factor, stroke_width=2)
            atom.set_fill(gold_color, opacity=1.0)
            atom.set_stroke(GOLD_E, width=2)
            # Scale and shift to center
            scaled_pos = (pos - center_offset) * scale_factor + np.array([0, -0.5, 0])
            atom.move_to(scaled_pos)
            atoms_final_2.add(atom)
        
        # Bonds for Part 2
        for i in range(lattice_size):
            for j in range(lattice_size):
                idx = i * lattice_size + j
                atom_pos = atoms_final_2[idx].get_center()
                if j < lattice_size - 1:
                    neighbor_pos = atoms_final_2[idx + 1].get_center()
                    bond = Line(atom_pos, neighbor_pos, color=WHITE, stroke_width=1.5)
                    bonds_final_2.add(bond)
                if i < lattice_size - 1:
                    neighbor_pos = atoms_final_2[idx + lattice_size].get_center()
                    bond = Line(atom_pos, neighbor_pos, color=WHITE, stroke_width=1.5)
                    bonds_final_2.add(bond)
        
        label_2 = Text("EAM", font_size=20, color=cloud_color)
        label_2.next_to(atoms_final_2, UP, buff=0.3)
        
        # Part 3 lattice (ML) - RIGHT
        atoms_final_3 = VGroup()
        bonds_final_3 = VGroup()
        for i, pos in enumerate(ml_positions):
            atom = Circle(radius=atom_radius * scale_factor, stroke_width=2)
            atom.set_fill(gold_color, opacity=1.0)
            atom.set_stroke(GOLD_E, width=2)
            # Scale and shift to right
            scaled_pos = (pos - center_offset) * scale_factor + np.array([3.5, -0.5, 0])
            atom.move_to(scaled_pos)
            atoms_final_3.add(atom)
        
        # Bonds for Part 3
        for i in range(lattice_size):
            for j in range(lattice_size):
                idx = i * lattice_size + j
                atom_pos = atoms_final_3[idx].get_center()
                if j < lattice_size - 1:
                    neighbor_pos = atoms_final_3[idx + 1].get_center()
                    bond = Line(atom_pos, neighbor_pos, color=WHITE, stroke_width=1.5)
                    bonds_final_3.add(bond)
                if i < lattice_size - 1:
                    neighbor_pos = atoms_final_3[idx + lattice_size].get_center()
                    bond = Line(atom_pos, neighbor_pos, color=WHITE, stroke_width=1.5)
                    bonds_final_3.add(bond)
        
        label_3 = Text("ML", font_size=20, color=RED)
        label_3.next_to(atoms_final_3, UP, buff=0.3)
        
        # Show all three lattices simultaneously
        self.play(
            FadeIn(bonds_final_1),
            FadeIn(atoms_final_1),
            FadeIn(label_1),
            FadeIn(bonds_final_2),
            FadeIn(atoms_final_2),
            FadeIn(label_2),
            FadeIn(bonds_final_3),
            FadeIn(atoms_final_3),
            FadeIn(label_3),
            run_time=2.0
        )
        self.wait(3.0)
        
        # Fade out everything
        self.play(
            FadeOut(title),
            FadeOut(comparison_label),
            FadeOut(bonds_final_1),
            FadeOut(atoms_final_1),
            FadeOut(label_1),
            FadeOut(bonds_final_2),
            FadeOut(atoms_final_2),
            FadeOut(label_2),
            FadeOut(bonds_final_3),
            FadeOut(atoms_final_3),
            FadeOut(label_3),
            run_time=1.5
        )
        self.wait(1.0)