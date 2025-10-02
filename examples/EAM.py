from manim import *
import numpy as np


class EAM(Scene):
    """Embedded Atom Method visualization with 5x5 gold lattice."""
    #Needs some work, vectors currently not updating. Could also flip the red arrows to better show repulsion.
    #Maybe we can actually start with all of the different arrows and combine them into one repulsion arrow due to the cloud.
    # Figure out how to get the arrows to update with the rest of the atom movement.
    
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
        # ALL 24 atoms will have forces with the missing atom
        neighbor_positions = lattice_positions  # Use all positions, not just nearest neighbors
        
        # Get the missing atom's current position
        missing_atom_pos = missing_atom.get_center()
        
        # Create force vectors pointing to/from the missing atom
        pair_potential_arrows = VGroup()
        electron_repulsion_arrows = VGroup()
        
        for neighbor_pos in neighbor_positions:
            # Direction from neighbor to missing atom position
            direction = missing_atom_pos - neighbor_pos
            distance = np.linalg.norm(direction)
            direction_norm = direction / distance
            
            # Scale arrow size by distance (closer = stronger force = longer arrow)
            distance_factor = np.clip(1.5 / (distance / spacing), 0.3, 1.0)
            
            # Pair potential force (attractive) - from each lattice atom toward missing atom
            pair_length = 0.4 * distance_factor
            pair_start = neighbor_pos + direction_norm * (atom_radius + 0.05)
            pair_end = neighbor_pos + direction_norm * (atom_radius + 0.05 + pair_length)
            pair_arrow = Arrow(
                start=pair_start,
                end=pair_end,
                color=pair_potential_color,
                buff=0,
                stroke_width=3 * distance_factor,
                max_tip_length_to_length_ratio=0.2
            )
            pair_potential_arrows.add(pair_arrow)
            
            # Electron cloud repulsion (repulsive) - from missing atom toward each lattice atom
            repulsion_direction = -direction_norm
            repulsion_length = 0.25 * distance_factor
            repulsion_start = missing_atom_pos + repulsion_direction * (cloud_radius * 0.3)
            repulsion_end = missing_atom_pos + repulsion_direction * (cloud_radius * 0.3 + repulsion_length)
            repulsion_arrow = Arrow(
                start=repulsion_start,
                end=repulsion_end,
                color=electron_repulsion_color,
                buff=0,
                stroke_width=2.5 * distance_factor,
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
        
        # Create updaters for force vectors to follow the atom and update directions
        def update_pair_arrows(mob):
            missing_pos = missing_atom.get_center()
            
            for i, neighbor_pos in enumerate(neighbor_positions):
                if i >= len(mob):
                    break
                    
                direction = missing_pos - neighbor_pos
                distance = np.linalg.norm(direction)
                if distance < 0.01:  # Avoid division by zero
                    continue
                direction_norm = direction / distance
                
                # Scale arrow size by distance (closer = stronger force)
                distance_factor = np.clip(1.5 / (distance / spacing), 0.3, 1.0)
                pair_length = 0.4 * distance_factor
                
                # Update attractive arrows (from lattice atoms to missing atom)
                new_start = neighbor_pos + direction_norm * (atom_radius + 0.05)
                new_end = neighbor_pos + direction_norm * (atom_radius + 0.05 + pair_length)
                new_arrow = Arrow(
                    start=new_start,
                    end=new_end,
                    color=pair_potential_color,
                    buff=0,
                    stroke_width=3 * distance_factor,
                    max_tip_length_to_length_ratio=0.2
                )
                mob[i].become(new_arrow)
        
        def update_repulsion_arrows(mob):
            missing_pos = missing_atom.get_center()
            
            for i, neighbor_pos in enumerate(neighbor_positions):
                if i >= len(mob):
                    break
                    
                direction = missing_pos - neighbor_pos
                distance = np.linalg.norm(direction)
                if distance < 0.01:  # Avoid division by zero
                    continue
                direction_norm = direction / distance
                
                # Scale arrow size by distance
                distance_factor = np.clip(1.5 / (distance / spacing), 0.3, 1.0)
                repulsion_length = 0.25 * distance_factor
                
                # Update repulsive arrows (from missing atom to lattice atoms)
                repulsion_direction = -direction_norm
                new_start = missing_pos + repulsion_direction * (cloud_radius * 0.3)
                new_end = missing_pos + repulsion_direction * (cloud_radius * 0.3 + repulsion_length)
                new_arrow = Arrow(
                    start=new_start,
                    end=new_end,
                    color=electron_repulsion_color,
                    buff=0,
                    stroke_width=2.5 * distance_factor,
                    max_tip_length_to_length_ratio=0.25
                )
                mob[i].become(new_arrow)
        
        # Add updaters to arrow groups
        pair_potential_arrows.add_updater(update_pair_arrows)
        electron_repulsion_arrows.add_updater(update_repulsion_arrows)
        
        # Animate atom moving into center with forces updating dynamically
        self.play(
            missing_atom.animate.move_to(center_pos),
            missing_cloud.animate.move_to(center_pos),
            run_time=2.5,
            rate_func=smooth
        )
        
        # Remove updaters and fade out forces
        pair_potential_arrows.remove_updater(update_pair_arrows)
        electron_repulsion_arrows.remove_updater(update_repulsion_arrows)
        self.play(
            FadeOut(pair_potential_arrows),
            FadeOut(electron_repulsion_arrows),
            FadeOut(force_label_pair),
            FadeOut(force_label_electron),
            run_time=0.5
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
