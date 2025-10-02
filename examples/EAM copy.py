from manim import *
import numpy as np


class EAM(Scene):
    """Embedded Atom Method visualization with 5x5 gold lattice."""
    
    def construct(self):
        # Lattice parameters
        lattice_size = 5
        spacing = 1.1
        atom_radius = 0.18
        cloud_radius = 0.52  # Reaches approximately to next atom
        
        # Colors
        gold_color = GOLD
        cloud_color = YELLOW
        pair_potential_color = GREEN  # Attractive
        electron_repulsion_color = RED  # Repulsive
        
        # Apply pseudo-3D isometric transformation
        def apply_perspective(pos):
            """Apply isometric projection for 3D look."""
            x, y = pos
            # Isometric projection
            iso_x = x - y * 0.5
            iso_y = (x + y) * 0.3
            return np.array([iso_x, iso_y, 0])
        
        # Create lattice positions
        lattice_positions = []
        center_pos = None
        center_idx = (lattice_size // 2, lattice_size // 2)
        
        for i in range(lattice_size):
            for j in range(lattice_size):
                x = (i - lattice_size // 2) * spacing
                y = (j - lattice_size // 2) * spacing
                pos_2d = np.array([x, y])
                pos = apply_perspective(pos_2d)
                
                if (i, j) == center_idx:
                    center_pos = pos
                else:
                    lattice_positions.append(pos)
        
        # Create atoms and electron clouds
        atoms = VGroup()
        electron_clouds = VGroup()
        
        for pos in lattice_positions:
            # Electron cloud (semi-transparent circle with gradient)
            cloud = Circle(radius=cloud_radius, stroke_width=0)
            cloud.set_fill(cloud_color, opacity=0.15)
            cloud.set_stroke(cloud_color, width=1, opacity=0.3)
            cloud.move_to(pos)
            electron_clouds.add(cloud)
            
            # Gold atom (circle with radial gradient effect)
            atom = Circle(radius=atom_radius, stroke_width=3)
            atom.set_fill(gold_color, opacity=1.0)
            atom.set_stroke(GOLD_E, width=3)
            atom.move_to(pos)
            atoms.add(atom)
        
        # Create the missing center atom (starts off-screen above)
        missing_cloud = Circle(radius=cloud_radius, stroke_width=0)
        missing_cloud.set_fill(cloud_color, opacity=0.15)
        missing_cloud.set_stroke(cloud_color, width=1, opacity=0.3)
        missing_cloud.move_to(center_pos + UP * 3.5)
        
        missing_atom = Circle(radius=atom_radius, stroke_width=3)
        missing_atom.set_fill(gold_color, opacity=1.0)
        missing_atom.set_stroke(GOLD_E, width=3)
        missing_atom.move_to(center_pos + UP * 3.5)
        
        # Title
        title = Text("Embedded Atom Method", font_size=36)
        title.to_edge(UP)
        
        # Add lattice
        self.add(title)
        self.play(
            LaggedStart(*[FadeIn(cloud, scale=0.5) for cloud in electron_clouds], lag_ratio=0.03),
            run_time=1.5
        )
        self.play(
            LaggedStart(*[GrowFromCenter(atom) for atom in atoms], lag_ratio=0.03),
            run_time=1.5
        )
        
        self.wait(0.5)
        
        # Highlight the vacancy
        vacancy_indicator = Circle(radius=cloud_radius * 1.2, color=WHITE, stroke_width=4)
        vacancy_indicator.move_to(center_pos)
        vacancy_label = Text("Vacancy", font_size=24, color=WHITE)
        vacancy_label.next_to(vacancy_indicator, DOWN, buff=0.3)
        
        self.play(Create(vacancy_indicator), FadeIn(vacancy_label, shift=UP * 0.2))
        self.wait(0.8)
        self.play(FadeOut(vacancy_indicator), FadeOut(vacancy_label))
        
        # Show the missing atom above
        self.play(
            FadeIn(missing_atom, shift=DOWN * 0.3),
            FadeIn(missing_cloud, shift=DOWN * 0.3)
        )
        self.wait(0.5)
        
        # Calculate force vectors before atom moves
        # Find nearest neighbors to center position (4 direct neighbors)
        neighbor_positions = []
        
        for i, pos in enumerate(lattice_positions):
            dist = np.linalg.norm(pos - center_pos)
            if 0.9 * spacing < dist < 1.3 * spacing:  # Nearest neighbors only
                neighbor_positions.append(pos)
        
        # Create force vectors from neighbors to center
        pair_potential_arrows = VGroup()
        electron_repulsion_arrows = VGroup()
        
        for neighbor_pos in neighbor_positions:
            # Direction from neighbor to center
            direction = center_pos - neighbor_pos
            direction_norm = direction / np.linalg.norm(direction)
            
            # Pair potential force (attractive) - from neighbor toward center
            pair_start = neighbor_pos + direction_norm * (atom_radius + 0.05)
            pair_end = neighbor_pos + direction_norm * (atom_radius + 0.45)
            pair_arrow = Arrow(
                start=pair_start,
                end=pair_end,
                color=pair_potential_color,
                buff=0,
                stroke_width=4,
                max_tip_length_to_length_ratio=0.25
            )
            pair_potential_arrows.add(pair_arrow)
            
            # Electron cloud repulsion (repulsive) - from center toward neighbor
            repulsion_direction = -direction_norm
            repulsion_start = center_pos + repulsion_direction * 0.05
            repulsion_end = center_pos + repulsion_direction * 0.5
            repulsion_arrow = Arrow(
                start=repulsion_start,
                end=repulsion_end,
                color=electron_repulsion_color,
                buff=0,
                stroke_width=4,
                max_tip_length_to_length_ratio=0.25
            )
            electron_repulsion_arrows.add(repulsion_arrow)
        
        # Show force vectors
        force_label_pair = Text("Pair Potential (Attractive)", font_size=22, color=pair_potential_color)
        force_label_pair.to_corner(UL, buff=0.5).shift(DOWN * 0.3)
        
        force_label_electron = Text("Electron Repulsion", font_size=22, color=electron_repulsion_color)
        force_label_electron.next_to(force_label_pair, DOWN, buff=0.25)
        
        self.play(Write(force_label_pair))
        self.play(
            LaggedStart(*[GrowArrow(arrow) for arrow in pair_potential_arrows], lag_ratio=0.15),
            run_time=1.5
        )
        self.wait(0.8)
        
        self.play(Write(force_label_electron))
        self.play(
            LaggedStart(*[GrowArrow(arrow) for arrow in electron_repulsion_arrows], lag_ratio=0.15),
            run_time=1.5
        )
        self.wait(1.5)
        
        # Animate atom moving into center
        self.play(
            missing_atom.animate.move_to(center_pos),
            missing_cloud.animate.move_to(center_pos),
            FadeOut(pair_potential_arrows),
            FadeOut(electron_repulsion_arrows),
            FadeOut(force_label_pair),
            FadeOut(force_label_electron),
            run_time=2.5,
            rate_func=smooth
        )
        
        # Add the atom to the complete lattice
        atoms.add(missing_atom)
        electron_clouds.add(missing_cloud)
        
        # Show complete lattice
        completion_text = Text("Complete Lattice", font_size=28, color=GREEN)
        completion_text.to_edge(DOWN)
        
        self.play(Write(completion_text))
        self.wait(1)
        
        # Pulse effect to show equilibrium
        self.play(
            atoms.animate.scale(1.1),
            electron_clouds.animate.scale(1.1),
            rate_func=there_and_back,
            run_time=1
        )
        
        self.wait(2)
